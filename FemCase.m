classdef FemCase < hgsetget
    properties
        mesh
        SElist  %% ???
        bc
        loads
        displacements
        reactions
        ncambio
        
        node_stress_array
    end
    properties (Dependent)
        nnodes
        nnel
        nodes_per_ele
        dofs_per_node
        dofs_per_ele
        n_node_dofs
        n_ele_dofs
        n_stress_components
        dim
    end
    
    methods (Static)
   
        function [global_coordinates, ncambio] = nlocal2global(coordinates,coordinates2)
            ncambio = zeros(size(coordinates2,1),1);
            tol = 1e-10;
            maxnodes = size(coordinates,1);
            count = 1;
            global_coordinates = coordinates;
            for i = 1:size(coordinates2,1)
                for j = 1:maxnodes
                    if norm(coordinates2(i,:) - coordinates(j,:)) < tol
                        ncambio(i) = j;
                        break
                    end
                end
                if j == maxnodes
                    ncambio(i) = maxnodes + count;
                    global_coordinates = [global_coordinates; coordinates2(i,:)];
                    count = count + 1;
                end
            end
        end
        function global_connections = global_connect(connections,ncambio)
            global_connections = zeros(size(connections));
            for i = 1:size(connections,1)
                for j = 1:size(connections,2)
                    global_connections(i,j) = ncambio(connections(i,j));
                end
            end
        end
        function obj =  case_from_SElist(SElist)
            SE1 = SElist(1);
            SE2 = SElist(2);
            [global_coordinates, ncambio_in] = FemCase.nlocal2global(SE1.get('coordinates'),SE2.get('coordinates'));
            global_connections = FemCase.global_connect(SE1.get('connections'),ncambio_in);
            mesh_fem = Mesh(global_coordinates,[SE1.get('connections');global_connections],SE1.get('material'));
            mesh_fem.set('ncambio',ncambio_in);
            obj = FemCase(mesh_fem,mesh_fem.ndofs(),mesh_fem.ndofs());
            SE1.set('ncambio',1:SE1.nnodes())
            SE2.set('ncambio',ncambio_in)
            obj.set('SElist',SElist)
        end
        function l_index1 = l_index(bsindex,index_out)
            l_index1 = zeros(length(index_out),1);
            for i = 1:length(index_out)
                l_index1(i) = find(bsindex == index_out(i));
            end
            
        end    
    end
    
    methods
        
        %% Initialize
        
        function femcase = FemCase(mesh_in,varargin)
            % femcase = FemCase(mesh_in)
            % femcase = FemCase(mesh_in,bc_in,loads_in)
            femcase.set('mesh',mesh_in);
            if isempty(varargin)
                BC_in = BC(femcase.nnodes,femcase.dofs_per_node, ...
                                femcase.nnel, femcase.dofs_per_ele);
                Loads_in = Loads(femcase.nnodes,femcase.dofs_per_node, ...
                            femcase.nnel, femcase.dofs_per_ele);
            else BC_in = varargin{1};
                Loads_in = varargin{2};
            end
            femcase.set('bc',BC_in);
            femcase.set('loads',Loads_in);
            Displacements_in = Displacements(femcase.nnodes, ...
                femcase.dofs_per_node, femcase.nnel, femcase.dofs_per_ele);
            femcase.set('displacements',Displacements_in);
            Reactions_in = Displacements(femcase.nnodes,femcase.dofs_per_node, ...
                            femcase.nnel, femcase.dofs_per_ele);
            femcase.set('reactions',Reactions_in);
        end
        %% Solve

        function solve_condensed(obj)
            % Takes into account superelements
            SElist = obj.get('SElist');
            SE1 = SElist(1);
            SE2 = SElist(2);
            
            index_in = SE1.nodesindex(SE1.nodes_in());
            index_out = SE1.nodesindex(SE1.nodes_out());
            index_in2 = SE2.nodesindex(SE2.nodes_in());
            index_out2 = SE2.nodesindex(SE2.nodes_out());
            
            stiffness = SE1.stiff();
            Kcc = stiffness(index_in,index_in);
            Krr = stiffness(index_out,index_out);
            Kcr = stiffness(index_in,index_out);
            Krc = Kcr';
            
            K_condensed = Krr- Krc*(Kcc\Kcr);
            
            nodos_unicos = unique([SE1.nodes_out() SE2.nodes_out()]);
            bsindex = SE1.nodesindex(nodos_unicos); %todos los grados de libertad del condensado.
            nodos_cond1 = FemCase.l_index(nodos_unicos,SE1.nodes_out());  %los nodos del 1er elemento respecto al condensado
            nodos_cond2 = FemCase.l_index(nodos_unicos,SE2.nodes_out());  % es decir, tiene indices mas chiquitos.
            
            l_index1 = SE1.nodesindex(nodos_cond1); %grados de libertd de esos indeces chicos.
            l_index2 = SE2.nodesindex(nodos_cond2);
            
            condensed_ndofs = length(bsindex);
            
            SmallStiff = sparse(condensed_ndofs,condensed_ndofs);
            
            SmallStiff(l_index1,l_index1) = SmallStiff(l_index1,l_index1) + K_condensed;
            SmallStiff(l_index2,l_index2) = SmallStiff(l_index2,l_index2) + K_condensed;
            
            condensed_free = obj.get('bc').get('nodelist');
            condensed_free = condensed_free(bsindex);
            
            %             condensed_loadsv = reshape(loads',[],1);
            condensed_loadsv = obj.get('loads').get('nodelist');
            condensed_loadsv = condensed_loadsv(bsindex);
            
            condensed_dis = sparse(condensed_ndofs,1);
            condensed_dis(condensed_free) = SmallStiff(condensed_free,condensed_free) \ condensed_loadsv(condensed_free);
            condensed_reactions = SmallStiff*condensed_dis;
            
            ndofs = obj.get('mesh').ndofs();
            displacements = zeros(ndofs,1);
            displacements(index_out) = condensed_dis(l_index1); %desplazamientos de caras del 1er elemento
            displacements(index_in) = -(Kcc\Kcr)*displacements(index_out);%desplazamientos interiores del 1er elemento
            displacements(index_out2) = condensed_dis(l_index2); %desplazamientos de caras del 2do elemento
            displacements(index_in2) = -(Kcc\Kcr)*displacements(index_out2);%desplazamientos interiores del 2do elemento
            
            reactions = zeros(size(displacements));
            reactions(bsindex) = condensed_reactions;
            
            obj.set('displacements',VectorField(displacements));
            obj.set('reactions',VectorField(reactions));
        end
        
        function solve(femcase)
            meshnow = femcase.get('mesh');
            stiffness = meshnow.stiff;
            D = zeros(meshnow.ndofs,1);
            free = femcase.get('bc').all_dofs;
            loadsv = femcase.get('loads').all_dofs;
            D(free) = stiffness(free,free)\loadsv(free);
            femcase.get('displacements').dof_list_in(D);
            R = stiffness*D;
            femcase.get('reactions').dof_list_in(R);
            % some sanity checks
            % BELONGS IN REACTIONS
            tol = 1e-15;
            for i = 1:3
                %sum(R(i:3:end))
                if sum(R(i:3:end))>tol
                    warning('sumR does not check')
                    sum(R(i:3:end))
                end
            end
        end
        %% Stress
        function stress_out = get.node_stress_array(femcase)
        % stress_out = get.node_stress_array(femcase)
        % Computes the stress at each node by element, and then averages
        % the values for stress_out array [nnodes x n_stress_components]
            if isempty(femcase.node_stress_array)
                dis = femcase.get('displacements').all_dofs;
                node_stress = cell(femcase.nnodes,1);
                for ele = 1:femcase.nnel
                    element = femcase.get('mesh').create_ele(ele,dis);
                    ele_stress = element.node_stress;
                    for n = 1:femcase.nodes_per_ele
                        node = element.node_id_list(n);
                        node_stress{node} = [node_stress{node}; ele_stress(n,:)];
                    end
                end
                stress_out = zeros(femcase.nnodes,size(node_stress{1},2));
                for node = 1:femcase.nnodes
                    if size(node_stress{node},1) == 1
                        stress_out(node,:) = node_stress{node};
                    else stress_out(node,:) = mean(node_stress{node});
                    end
                end
                femcase.node_stress_array = stress_out;
            else femcase.node_stress_array;
            end
        end
            
        %% Plot
        
        function plot_displacement(femcase,coord)
            scale = 1000;
            hold on
            mesh_aux = femcase.get('mesh');
            node_func = femcase.get('displacements').node_function;
            mesh_aux.plot_node_function(node_func(:,coord))
            hold off
        end
        
        %% Setters & Getters
        function n = get.nnodes(femcase)
            n = femcase.get('mesh').nnodes;
        end
        function n = get.nnel(femcase)
            n = femcase.get('mesh').nnel;
        end
        function n = get.nodes_per_ele(femcase)
            n = femcase.get('mesh').nodes_per_ele;
        end
        function n = get.dofs_per_node(femcase)
            n = femcase.get('mesh').dofs_per_node;
        end
        function n = get.dofs_per_ele(femcase)
            n = femcase.get('mesh').dofs_per_ele;
        end
        function n = get.n_node_dofs(femcase)
            n = femcase.get('mesh').n_node_dofs;
        end
        function n = get.n_ele_dofs(femcase)
            n = femcase.get('mesh').n_ele_dofs;
        end

        function n = get.n_stress_components(femcase)
            n = femcase.get('mesh').create_ele(1).n_stress_components;
        end
        function n = get.dim(femcase)
            n = femcase.get('mesh').create_ele(1).dim;
        end
    end
    

    
    
end
