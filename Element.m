classdef (Abstract) Element < hgsetget
    properties
        id
        nodes
        material
        parent
        dofs_per_ele
        dofs_per_node
        dof_dis             % Solutions in dof form
    end
    properties (Abstract, Constant)
        nnodes              % Number of nodes in the element
        node_connectivity   % Rows contain pairs of connected nodes
        face_connectivity   % Nodes of each face forming a closed polygon
        node_local_coords   % Local coordinates of each node for the N(xi,eta,mu)
    end
    properties (Dependent)
        dim
        n_nodes
        node_id_list
        n_nodes_per_face
        n_faces
        n_node_dofs
        n_ele_dofs
        dof_list
        node_dof_list
        ele_dof_list
        n_dofs              % Total number of dofs in the element
    end
    methods (Abstract)
        N_out = N(element,local_coords)
        dN_out = dN_dxi(element,local_coords)
        dN_out = dN_deta(element,local_coords)
        dN_out = dN_dmu(element,local_coords)
        dN_out = DN(element,local_coords)
    end
        
    methods
        function obj = Element(id_in,dofs_per_node_in,dofs_per_ele_in,nodes, ...
                                        material)
            set(obj,'id',id_in);
            set(obj,'nodes',nodes);
            set(obj,'material',material);
            set(obj,'dofs_per_node',dofs_per_node_in);
            set(obj,'dofs_per_ele',dofs_per_ele_in);
        end

        function facecoord = elementface(obj,coordnum,value)
            coordinat = obj.coordinates;
            AL = Face.Al();
            count = 1;
            for i = 1:3
                for j = [-1 1];
                    if and(i==coordnum,j==value)
                        facecoord = coordinat(AL(count,:),:);
                    end
                    count = count + 1;
                end
            end
        end
        function facenodesid = facenodesid(obj,coordnum,value)
            nodesids = obj.node_id_list();
            AL = Face.Al();
            count = 1;
            for i = 1:3
                for j = [-1 1];
                    if and(i==coordnum,j==value)
                        facenodesid = nodesids(AL(count,:));
                        break
                    end
                    count = count + 1;
                end
            end
        end
        function coordinat = coordinates(element)
            nodelist = element.get('nodes');
            coordinat = zeros(numel(nodelist),element.dim);
            for i = 1:numel(nodelist)
                nodei = nodelist(i);
                coordinat(i,:) = nodei.get('coordinates');
            end
        end

        function jac_out = jacobian(element,local_coords)
            % jac_out [3x3] = jacobian(obj,xi,eta,mu)
            for i = local_coords
                if or(i>1,i<-1)
                    error('wrong args for xi,eta, or mu')
                end
            end
            jac_out = element.DN(local_coords)*element.coordinates();
            if det(jac_out)<0
                error('det(Jacobian) < 0')
            end
        end
        function face_array = local_faces_from_nodes_id(element,nodes_in)
            % face_array = faces_from_nodes_id(element,nodes_in)
            % Returns that faces that the nodes_in generate in the element
            % in polygon form. Takes a subset from face_connectivity.
            local_ids = element.node_id_to_local(nodes_in);
            face_array = element.faces_from_local_nodes(local_ids);
        end
        function face_array = global_faces_from_nodes_id(element,nodes_in)
            % face_array = global_faces_from_nodes_id(element,nodes_in)
            % Returns that faces that the nodes_in generate in the element
            % in polygon form in GLOBAL ids Takes a subset from 
            % face_connectivity and then transforms it to global ids.
            local_face_array = element.local_faces_from_nodes_id(nodes_in);
            face_array = element.node_id_to_global(local_face_array);
        end
        function face_array = faces_from_local_nodes(element,nodes_in)
            % face_array = faces_from_local_nodes(element,nodes)
            % Returns that faces that the nodes_in generate in the element
            % in polygon form. Takes a subset from face_connectivity.
            % nodes_in range from 1 to elemenet.nnodes
            require(max(nodes_in)<=element.nnodes,'Nodes List not in Local Ids')
            face_array = [];
            for f = 1:element.n_faces % Check in every face if all it's nodes are present
                n = 0;  % Node number I'm trying in the face (ranges from 1 to n_nodes_per_face
                node_here = true; % Checks if the node is in the face
                while (node_here) && n < element.n_nodes_per_face
                    n = n + 1;
                    node_here = any(element.face_connectivity(f,n)==nodes_in);
                end
                if node_here && n == element.n_nodes_per_face % All nodes where in the list!
                    face_array = [face_array; element.face_connectivity(f,:)];  % Add it
                end
            end
        end
        function global_ids = node_id_to_global(element,nodes_in)
            % global_ids = node_id_to_global(element,nodes_in)
            % nodes_in [m x n] -> global_ids [m x n]
            % Takes a list of local ids and returns their global
            % counterparts. Inverse of node_id_to_local
            % Every node in nodes_in should be < element.nnodes
            func = @(node) iif(node <= element.nnodes, ...
                                element.node_id_list(node), NaN);
            global_ids = arrayfun(func,nodes_in);
        end
        function local_ids = node_id_to_local(element,nodes)
            % local_ids = node_ids_to_local(element,nodes)
            % Takes a list of node ids, compares to the elements own node
            % ids and returns in local numbers the ones that it
            % contains. Might return an empty list if none are found.
            local_ids = [];
            for n = 1:length(nodes)
                local_ids = [local_ids find(nodes(n)==element.node_id_list)];
            end
        end
        function cn_out = connected_nodes(element,nodenum)
            % cn_out = connected_nodes(element,nodenum)
            % returns a 1x? vector with the nodes nodenum is connected to.
            require(nodenum<=element.nnodes,'Input node bigger than number of nodes');
            NC = element.node_connectivity;
            cn_out = [];
            for n = 1:element.nnodes
                if any(NC(n,:) == nodenum)
                    cn_out = [cn_out NC(n,NC(n,:)~=nodenum)];
                end
            end
        end
        
        
        %% DOFS
        
        function dof_list = get.dof_list(element)
            % dof_list = dof_list(mesh,ele_num)
            % Returns a dof_list [nodes_per_element*n_dof_per_node +
            % n_dofs_per_element, 1] for the element number ele_num
            % Element DOFs are always after the node dofs.
            dof_list = [element.node_dof_list; element.ele_dof_list];
        end
        function dof_list = get.node_dof_list(element)
            % dof_list = node_dof_list(mesh,ele_num)
            % Returns a list with all the element - node_dofs
            dof_list = zeros(element.nnodes*element.dofs_per_node,1);
            for n = 1:element.nnodes;
                node = element.node_id_list(n);
                dof_list(index_range(element.dofs_per_node,n)) =  ...
                                        index_range(element.dofs_per_node,node);
            end
        end
        function dof_list = get.ele_dof_list(element)
            % dof_list = dof_list(mesh,ele_num)
            % Returns a list with all the element - element_dofs
            dof_list = index_range(element.n_ele_dofs,element.get('id'));
        end

        %% Setters and Getters
        
        function [xi, eta, mu] = components(element,local_coords)
            if length(local_coords) == 3
                xi = local_coords(1);
                eta = local_coords(2);
                mu = local_coords(3);
            elseif length(local_coords) == 2
                xi = local_coords(1);
                eta = local_coords(2);
                if element.dim == 3
                    mu = 0;
                else mu = NaN;
                end
%             else error('Dimension missmatch local_coords and element dim')
            end
                
        end
        function nodeids = get.node_id_list(obj)
        	nodelist = obj.get('nodes');
            nodeids = zeros(size(nodelist));
            for i = 1:numel(nodelist)
                node = nodelist(i);
                nodeids(i) = node.get('id');
            end
        end
        function dim_out = get.dim(element)
            node_list = element.get('nodes');
            dim_out = length(node_list(1).get('coordinates'));
        end
        function dim_out = get.n_nodes(element)
            dim_out = length(element.get('nodes'));
        end
        function nnodes = get.n_nodes_per_face(element)
            nnodes = size(element.get('face_connectivity'),2);
        end
        function nnodes = get.n_faces(element)
            nnodes = size(element.get('face_connectivity'),1);
        end
        function dofs_n = get.n_node_dofs(element)
            dofs_n = element.nnodes*element.dofs_per_node;
        end
        function dofs_n = get.n_ele_dofs(element)
            dofs_n = element.dofs_per_ele;
        end
        function dofs_n = get.n_dofs(element)
            dofs_n = element.n_ele_dofs + element.n_node_dofs;
        end
%         function dof_n = get.last_node_dof(element)    % Last dof belonging to a node
%             dof_n = element.dofs_per_node*element.nnodes;
%         end
        
        function set.id(element,id_in)
            require(isscalar(id_in) && mod(id_in,1)==0,'Wrong element Id');
            element.id = id_in;
        end
        function set.nodes(element,nodes_in)
            require(length(nodes_in)==element.nnodes,'Wrong amount of Nodes');
            require(isa(nodes_in(1),'Node'),'Node List doesn''t contain nodes');
            element.nodes = nodes_in;
        end
        function set.material(element,material_in)
            require(isa(material_in,'Material'),'Wrong Input Type for Material');
            element.material = material_in;
        end
        function set.dof_dis(element,dis_in)
        	require(length(dis_in)==element.n_dofs,'Wrong Displacement Size');
            element.dof_dis = dis_in;
        end
    end
 

end
