classdef FemCase < hgsetget
    properties
        mesh
        SElist
        bc
        loads
        displacements
        reactions
        ncambio
        StressArray
        StrainArray
        StressArrayEle 
        StrainArrayEle
        MaxStress
        MaxStressInd
        MaxStressPoint
        MinStress
        MinStressInd
        MinStressPoint
        MaxStressEle
        MaxStressIndEle
        MaxStressPointEle
        MinStressEle
        MinStressIndEle
        MinStressPointEle
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
            obj = FemCase(mesh_fem,mesh_fem.sdof(),mesh_fem.sdof());
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
        
        function femcase = FemCase(mesh_in,bc_in,loads_in)
            femcase.set('mesh',mesh_in);
            if isa(bc_in,'BC')
                BC_in = bc_in;
            elseif isempty(bc_in)
                BC_in = BC(femcase.nnodes,femcase.n_node_dofs, ...
                            femcase.nelements, femcase.n_element_dofs);
            end
            femcase.set('bc',BC_in);
            if isa(loads_in,'Loads')
                Loads_in = loads_in;
            elseif isempty(bc_in)
                Loads_in = Loads(femcase.nnodes,femcase.n_node_dofs, ...
                            femcase.nelements, femcase.n_elemens_dofs);
            end
            femcase.set('loads',Loads_in);
            Displacements_in = Displacements(femcase.nnodes,femcase.n_node_dofs, ...
                            femcase.nelements, femcase.n_element_dofs);
            femcase.set('displacements',Displacements_in);
            Reactions_in = Displacements(femcase.nnodes,femcase.n_node_dofs, ...
                            femcase.nelements, femcase.n_element_dofs);
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
            
            condensed_sdof = length(bsindex);
            
            SmallStiff = sparse(condensed_sdof,condensed_sdof);
            
            SmallStiff(l_index1,l_index1) = SmallStiff(l_index1,l_index1) + K_condensed;
            SmallStiff(l_index2,l_index2) = SmallStiff(l_index2,l_index2) + K_condensed;
            
            condensed_free = obj.get('bc').get('nodelist');
            condensed_free = condensed_free(bsindex);
            
            %             condensed_loadsv = reshape(loads',[],1);
            condensed_loadsv = obj.get('loads').get('nodelist');
            condensed_loadsv = condensed_loadsv(bsindex);
            
            condensed_dis = sparse(condensed_sdof,1);
            condensed_dis(condensed_free) = SmallStiff(condensed_free,condensed_free) \ condensed_loadsv(condensed_free);
            condensed_reactions = SmallStiff*condensed_dis;
            
            sdof = obj.get('mesh').sdof();
            displacements = zeros(sdof,1);
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
            stiffness = meshnow.stiff();
            D = zeros(meshnow.sdof(),1);
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
        function tensors(femcase)
            dis_aux = femcase.get('displacements').all_dofs;
            StressA = zeros(femcase.nelements,8,6);
            StrainA = StressA;
            StressA_Ele = zeros(femcase.nelements,6);
            StrainA_Ele = StressA_Ele;
            gaussn = 2;
            gaussp = lgwt(gaussn,-1,1);
            for ele = 1:femcase.nelements
                new_element = femcase.get('mesh').element_create(ele);
                count = 1;
                for i = 1:gaussn
                    xi = gaussp(i);
                    for j = 1:gaussn
                        eta = gaussp(j);
                        for k = 1:gaussn
                            mu = gaussp(k);
                            aux = new_element.B(xi,eta,mu)*dis_aux(new_element.ldofsid());
                            C = new_element.C;
                            StrainA(ele,count,:) = aux;
                            StressA(ele,count,:) = C*aux;
                            count = count + 1;
                            StrainA_Ele(ele,:) = StrainA_Ele(ele,:) + aux'/8;
                            StressA_Ele(ele,:) = StressA_Ele(ele,:) + (C*aux)'/8;
                        end
                    end
                end
            end
            femcase.set('StrainArray',StrainA);
            femcase.set('StressArray',StressA);
            femcase.set('StrainArrayEle',StrainA_Ele);
            femcase.set('StressArrayEle',StressA_Ele);
        end
        function ExtremeStress(obj,opts)
            coordinates = obj.get('mesh').get('coordinates');
            MaxStress = zeros(3,1);
            MaxStressInd = zeros(3,1);
            MaxStressPoint = zeros(3,3);
            for i = 1:3
                switch opts
                    case 1
                        [val, ind] = max(StressArray(:,:,i));
                        [val, ind2] = max(val);
                    case 2
                        [val, ind] = min(StressArray(:,:,i));
                        [val, ind2] = min(val);
                end     
                MaxStress(i) = val;
                MaxStressInd(i) = ind(ind2);
                MaxStressPoint(i,:) = coordinates(connections(MaxStressInd(i),1),:); %Devuelve la coordenada del nodo uno de ese elemento
            end
            switch opts
                case 1
                    obj.set('MaxStress',MaxStress)
                    obj.set('MaxStressInd',MaxStressInd)
                    obj.set('MaxStressPoint',MaxStressPoint)
                case 2
                    obj.set('MinStress',MaxStressEle)
                    obj.set('MinStressInd',MaxStressInd)
                    obj.set('MinStressPoint',MaxStressPoint)
            end
        end
        function ExtremeStressElement(obj,opts)
            StressArrayEle = obj.get('StressArrayEle');
            StrainArrayEle = obj.get('StrainArrayEle');
            coordinates = obj.get('mesh').get('coordinates');
            connections = obj.get('mesh').get('connections');
            MaxStressEle = zeros(3,1);
            MaxStressIndEle = zeros(3,1);
            MaxStressPointEle = zeros(3,3);
            for i = 1:3
                switch opts
                    case 1
                        [val ind] = max(StressArrayEle(:,i));
                    case 2
                        [val ind] = min(StressArrayEle(:,i));
                end   
                MaxStressEle(i) = val;
                MaxStressIndEle(i) = ind;
                MaxStressPointEle(i,:) = coordinates(connections(ind,1),:); %Devuelve la coordenada del nodo uno de ese elemento
            end
            switch opts
                case 1
                    obj.set('MaxStressEle',MaxStressEle)
                    obj.set('MaxStressIndEle',MaxStressIndEle)
                    obj.set('MaxStressPointEle',MaxStressPointEle)
                case 2
                    obj.set('MinStressEle',MaxStressEle)
                    obj.set('MinStressIndEle',MaxStressIndEle)
                    obj.set('MinStressPointEle',MaxStressPointEle)
            end
        end
        
        function getMaxStressEle(obj)
            obj.tensors();
            obj.ExtremeStressElement(1)
        end
            
        %% Plot
        
        function plot(obj,elelist)
            scale = 1000;
            hold on
            mesh_aux = obj.get('mesh');
            mesh_aux = Mesh(mesh_aux.get('element_type'),mesh_aux.get('coordinates'), ...
                    mesh_aux.get('connections'),mesh_aux.get('material'));
            dis = obj.get('displacements');
            Dis = dis.node_function;
            mesh_aux.plot('g');
            bc_aux = obj.get('bc');
            bc_aux = VectorField(~bc_aux.get('nodelist'));
            bc_aux.plotVF(mesh_aux.get('coordinates'),'k')
            if ~isempty(dis)
                dis.plotVF(mesh_aux.get('coordinates'),'b')
            end
            obj.get('loads').plotVF(mesh_aux.get('coordinates'),'r');
            hold off
        end
        
        %% Helpers
        
        function n = n_node_dofs(femcase)
            n = femcase.get('mesh').get('n_node_dofs');
        end
        function n = n_element_dofs(femcase)
            n = femcase.get('mesh').get('n_element_dofs');
        end
        function n = nnodes(femcase)
            n = femcase.get('mesh').nnodes;
        end
        function n = nelements(femcase)
            n = femcase.get('mesh').nnel;
        end
        
    end
    
    
end
