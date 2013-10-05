classdef Mesh < hgsetget
    properties
        element_type
        coordinates
        connections
        material
        dofs_per_node
        dofs_per_ele
        outer_face
        polygons            % Array(:,:,1:3) = [X;Y;Z] for patch. Stored for Speed
        parent
        ncambio
    end
    
    properties (Dependent)
        nnodes
        nnel
        nodes_per_ele
        n_node_dofs         % Total dofs belonging to nodes
        n_ele_dofs          % Total dofs belonging to element but not nodes
        ndofs               % Total number of dofs in the mesh
        coords
        connect
        dim
        outer_nodes
        inner_nodes
        last_node_dof
    end
    
    methods

        function mesh = Mesh(type,coordinates,connections,material)
            set(mesh,'coordinates',coordinates);
            set(mesh,'connections',connections);
            set(mesh,'material',material);
            set(mesh,'element_type',type);
            ele = mesh.create_ele(1);
            set(mesh,'dofs_per_node',ele.get('dofs_per_node'));       % Dofs per node
            set(mesh,'dofs_per_ele',ele.get('dofs_per_ele'));  
        end

        %% Coordinates Methods
        
        function value = mincoordval(mesh,coordnum)
            coordinatesl = mesh.get('coordinates');
            [~, ind] = min(abs(coordinatesl(:,coordnum)));
            value = coordinatesl(ind,coordnum);
        end
        function value = maxcoordval(mesh,coordnum)
            coordinatesl = mesh.get('coordinates');
            [~, ind] = max(abs(coordinatesl(:,coordnum)));
            value = coordinatesl(ind,coordnum);
        end
        function shift(mesh,dis_vector)
            coords = mesh.get('coordinates');
            coords2 = zeros(size(coords));
            for n = 1:mesh.nnodes()
                coords2(n,:) = coords(n,:) + dis_vector;
            end
            mesh.set('coordinates',coords2)
        end
        
        %% Face Methods
        
        function add_outer_face(mesh,newface)
            mesh.set('outer_face',[mesh.get('outer_faces');newface])
        end

        %% Node Methods
        function nodesin = nodes_in(mesh)
            nodesin = 1:size(mesh.get('coordinates'),1);      
            nodesin(mesh.get('nodesout')) = [];   
            ncambio_aux = mesh.get('ncambio');
            nodesout = reshape(ncambio_aux(nodesin),1,[]); % WHAT IS HERE?
        end
        function nodesout = nodes_out(mesh)
            ncambio_aux = mesh.get('ncambio');
            nodesout = reshape(ncambio_aux(mesh.get('nodesout')),1,[]);      
        end 
        function ind = nodesindex(mesh,nodelist)
            ndofpernode = mesh.dofs_per_node();
            ind = zeros(length(nodelist),1);
            for i = 1:length(nodelist)
                range = (ndofpernode*(i-1)+1):(ndofpernode*i);
                ind(range) = (ndofpernode*(nodelist(i)-1)+1):(ndofpernode*nodelist(i));
            end
        end
        function ele_list = elementsofnode(mesh,node_in)
            % elementlist = elementsofnode(mesh,nodenum)
            % Returns a list of element numbers that contain node_in
            ele_list = [];
            for ele = 1:mesh.nnel
                if any(node_in == mesh.connect(ele,:))
                    ele_list = [ele_list ele];
                end
            end
        end
        function random_inner(mesh,delta)
            % randominner(mesh,delta)
            % !!! state altering.
            % Changes the position of a random inner node by delta
            innodes = mesh.inner_nodes;
            if ~isempty(innodes)
                coord = mesh.get('coordinates');
                coord(innodes(1),1) = coord(innodes(1),1) + delta;
                mesh.set('coordinates',coord);
            end
        end
        function nodenum = findnode(mesh,x0)
            % nodenum = findnode(mesh,x0)
            % Finds a node close to x0 according to a tolerance.
            require(length(x0)==mesh,dim,'Wrong coordinates size');
            nodenum = [];
            tol = 5e-3;
            for n = 1:mesh.nnodes()
                if norm(mesh.coords(n,:)-x0)<tol;
                    nodenum = n;
                end
            end
        end
        %% Element
        function elementcoord_out = elementcoord(mesh,ele)
            % elementcoord_out = elementcoord(mesh,ele)
            % Find the center of gravity of the element and return it as 
            % a coordinate
            % SHOULD BE IN ELEMENTS
            new_ele = mesh.create_ele(ele);
            elementcoord_out = median(new_ele.coordinates());
        end
        
        function dof_list = ele_dof_list(mesh,ele_num)
            % dof_list = ele_dof_list(mesh,ele_num)
            % Returns a dof_list [nodes_per_element*n_dof_per_node +
            % n_dofs_per_element, 1] for the element number ele_num
            % Element DOFs are always after the node dofs.
            dof_list = [mesh.ele_node_dof_list(ele_num); mesh.ele_ele_dof_list];
        end
        function dof_list = ele_node_dof_list(mesh,ele_num)
            % dof_list = ele_node_dof_list(mesh,ele_num)
            % Returns a list with all the element - node_dofs
            element = mesh.create_ele(ele_num);
            dof_list = zeros(element.nnodes*element.n_node_dofs,1);
            for n = element.nnodes;
                node = element.node_id_list(n);
                dof_list(index_range(element.n_node_dofs,n)) =  ...
                                        index_range(element.n_node_dofs,node);
            end
        end
        function dof_list = ele_ele_dof_list(mesh,ele_num)
            % dof_list = ele_dof_list(mesh,ele_num)
            % Returns a list with all the element - element_dofs
            dof_list = index_range(mesh.n_ele_dofs,ele_num) + ...
                                    mesh.last_node_dof;
        end
        %% Assembly Methods
        function element = create_ele(mesh,ele_num,varargin)
            % element = create_ele(mesh,ele_num)
            % element = create_ele(mesh,ele_num,displacements)
            % CREATES THE APPROPIATE ELEMENT FROM ID: ele_num
            % If the displacements where found, then dis_in is stored
            connect = mesh.get('connections');
            coords = mesh.get('coordinates');
            for i = 1:mesh.nodes_per_ele;
                nodos(i) = Node(connect(ele_num,i),coords(connect(ele_num,i),:),mesh);
            end
            % BIG PROBLEM HERE -> How to initialize by element type?
            switch mesh.element_type
                case 'Mech_H8'
                    element = Mech_H8(ele_num,nodos,mesh.get('material'));   % HARDCODED
                case 'Mech_Q4'
                    element = Mech_Q4(ele_num,nodos,mesh.get('material'));
                otherwise
                    error(strcat('Element Type: ',mesh.element_type,' not found'));
            end
            if ~isempty(varargin)
                dis_aux = varargin{1};
                if length(dis_aux) == mesh.ndofs
                    dis_in = dis_aux(element.dof_list);
                else dis_in = dis_aux;
                end
                element.set('dof_dis',dis_in);
            end
        end
        function stiffness = stiff(mesh)
                nnel = mesh.nnel();
                stiffness = sparse(mesh.ndofs(),mesh.ndofs());
                gaussn = 2;
                for ele = 1:nnel
                    element = mesh.create_ele(ele);
                    kele = element.K(gaussn);
                    ind = element.ldofsid();
                    stiffness(ind,ind) = stiffness(ind,ind) + kele;
                end
        end
        
        
        %% Plot
        
        function [X,Y,Z] = out_polygons(mesh)
            % poly = out_polygons(mesh)
            % Returns a list of polygons for the MATLAB patch function
            poly_aux = mesh.get('polygons');
            X = poly_aux(:,:,1); Y = poly_aux(:,:,2); 
            if mesh.dim == 2
                Z = zeros(size(X));
            else Z = poly_aux(:,:,3);
            end
        end
        
        function plot(mesh,colour)
            hold on
            [X,Y,Z] = mesh.out_polygons;
            patch(X,Y,Z,colour)
            hold off
        end
        
        function plot_node_function(mesh,node_list)
            % plot_node_function(mesh,node_list)
            % Takes a list of values associated to nodes and plots them on
            % the mesh
            require(numel(node_list) == mesh.nnodes,'Wrong Size node_list')
            node_list = reshape(node_list,[],1);
            hold on
            [X,Y,Z] = mesh.out_polygons;
            func = @(node) node_list(node);
            C = mesh.node_fun_to_poly(func);
            require(all(size(X)==size(C)),'Error with node_fun_to_poly');
            patch(X,Y,Z,C)
            hold off
        end
        
        function poly_out = node_fun_to_poly(mesh,func)
            % poly_out = node_fun_to_poly(mesh,func)
            % Takes a node_function, that accepts node indexes and puts it
            % in outer polygon shape for patch            
            out_nodes = mesh.outer_nodes;
            poly_out = [];
            for ele = 1:mesh.nnel   % Loop every element
                ele_instance = mesh.create_ele(ele);
                % Get the outer faces of the element
                faces = ele_instance.global_faces_from_nodes_id(out_nodes)';
                for f = 1:size(faces,2) % Loop through the surfaces
                    % Store the coordinates of each outer node in order
                    poly_out = [poly_out func(faces(:,f))];
                end
            end
        end
        
        %% Dependet but Complex Property Methods
        function out_face = get.outer_face(mesh)
            % out_face = get.outer_face(mesh)
            % Returns a Face Object with  all the nodes that lie in the Boundary
            % It finds it by finding which elements are connected to less
            % elements.
            if isempty(mesh.outer_face)
                if mesh.dim == 2
                    out_face = Face(1:mesh.nnodes,mesh);
                    mesh.set('outer_face',out_face);
                else
                    nodelist = [];
                    %                 par_plot = [];
                    for nodenum = 1:mesh.nnodes()   % Loop through all the nodes
                        elelist = mesh.elementsofnode(nodenum);
                        % Finds the elements the node belongs too
                        element = mesh.create_ele(elelist(1));
                        if length(elelist) < element.n_nodes    % If
                            nodelist = [nodelist nodenum];
                            %                     for i = 1:length(elelist)
                            %                     	allnodes = mesh.connect(elelist(i),:);  %should be [n_nodesx1]
                            %                         local_node = find(nodenum == allnodes);
                            %                         other_nodes = allnodes(element.connected_nodes(local_node));
                            %                         for j = 1:length(other_nodes)
                            %                             if any(other_nodes(j)==nodelist)
                            %                                 par_plot = [par_plot; nodenum other_nodes(j)];
                            %                             end
                            %                         end
                            %                     end
                        end
                    end
                    %             big_elelist = unique(big_elelist);
                    out_face = Face(nodelist,mesh);
                    mesh.set('outer_face',out_face);
                    %             outer_faces.set('par_plot',par_plot);
                end
            else out_face = mesh.outer_face;
            end
        end
        
        function poly_out = get.polygons(mesh)
            if isempty(mesh.polygons)
                X = mesh.node_fun_to_poly(@(nodes) mesh.coords(nodes,1));
                poly_out = zeros([size(X) mesh.dim]);
                poly_out(:,:,1) = X; 
                for d = 2:mesh.dim
                    func = @(nodes) mesh.coords(nodes,d);
                    poly_out(:,:,d) = mesh.node_fun_to_poly(func); 
                end
                mesh.set('polygons',poly_out);
            else poly_out = mesh.polygons;
            end
        end
        

        %% Property Methods
        function in_nodes = get.inner_nodes(mesh)
            % in_nodes = get.inner_nodes(mesh)
            % Finds a list of inner node by comparing to the outer nodes
            % Computationally expensive, maybe stored for efficiency (???)
            in_nodes = [];
            count = 1;
            for n = 1:mesh.nnodes()
                if ~any(mesh.outer_nodes == n)
                    in_nodes(count) = n;
                    count = count + 1;
                end
            end
        end
        function out_nodes = get.outer_nodes(mesh)
            out_nodes = mesh.get('outer_face').get('node_list');
        end
        function out = get.nnodes(mesh)          %total number of nodes
            out = size(mesh.coords,1);
        end
        function out = get.nnel(mesh)            %total number of elements
            out = size(mesh.connect,1);
        end
        function out = get.nodes_per_ele(mesh)
            out = size(mesh.connect,2);
        end
        function out = get.n_node_dofs(mesh)    %number of dofs per node
            out = mesh.get('dofs_per_node')*mesh.nnodes;
        end
        function out = get.n_ele_dofs(mesh)    %number of dofs per node
            out = mesh.get('dofs_per_ele')*mesh.nnel;
        end
        function out = get.ndofs(mesh)            %total system dofs
            out = mesh.nnodes()*mesh.dofs_per_node();
        end
        function out = get.dim(mesh)    %number of dofs per node
            out = size(mesh.coords,2);
        end
        function coords_out = get.coords(mesh)
            coords_out = mesh.get('coordinates');
        end
        function coords_out = get.connect(mesh)
            coords_out = mesh.get('connections');
        end
        function dof_n = get.last_node_dof(mesh)    % Last dof belonging to a node
            dof_n = mesh.n_node_dofs*mesh.nnodes;
        end
    end  
    
     methods (Static)
        function L = NXorder()
            L = [6 2 5 1 7 3 8 4];
        end
        function sca = tometer()
            sca = 1e3;
        end
        function mesh = mesh_import(coordstr,connectstr,mate)
            connections_in = load(connectstr);
            connections_in = connections_in(:,2:end);
            coordinates_in = load(coordstr);
            connections2 = zeros(size(connections_in));
            for i = 1:size(coordinates_in,1)
                nodenum = coordinates_in(i,1);
                connections2(connections_in==nodenum) = i;
            end
            connections_in = connections2;
            coordinates_in = coordinates_in(:,2:4);
            l = Mesh.NXorder();
            aux = zeros(size(connections_in));
            for i = 1:length(l)
                aux(:,i) = connections_in(:,l(i));
            end
            connections_in = aux;
            coordinates_in = coordinates_in/1e4;
            coordinates_in = coordinates_in/Mesh.tometer();
            mesh = Mesh(coordinates_in,connections_in,mate);
        end
        function mesh = mesh_import2(coordstr,connectstr,mate)
            connections_in = load(connectstr);
            coordinates_in = load(coordstr);
            connections2 = zeros(size(connections_in));
            for i = 1:size(coordinates_in,1)
                nodenum = coordinates_in(i,1);
                connections2(connections_in==nodenum) = i;
            end
            
            coordinates_in = coordinates_in/1e4;
            coordinates_in = coordinates_in/Mesh.tometer();
            mesh = Mesh(coordinates_in,connections_in,mate);
        end
        function mesh2D = meshgen2D(type,A,M,E,nu,rho)
            % mesh2D = meshgen2D(A,M,E,nu,rho)
            % Generates a 2D rectangular mesh
            [coords_aux, connect_aux] = Mesh.rectangle_mesh([A 1],[M 1]);
            nnodes = size(coords_aux,1);
            material_in = Material(1,'Iso',E,nu,rho);
            mesh2D = Mesh(type,coords_aux(1:nnodes/2,1:2), ...
                            connect_aux(:,1:4),material_in);            
        end
        function mesh3D = meshgen(type,A,M,E,nu,rho)
            % mesh3D = meshgen(A,M,E,nu,rho)
            % Generates a 3D rectangular mesh
            [coordinates_in, connections_in] = Mesh.rectangle_mesh(A,M);
            mat1 = Material(1,'Iso',E,nu,rho);
            mesh3D = Mesh(type,coordinates_in,connections_in,mat1);
        end
        function [coordinates_out, connections_out] = rectangle_mesh(A,M)
            a = A(1);b = A(2); c = A(3);
            m = M(1);n = M(2); p = M(3);
            x = a/m;y = b/n;z = c/p;
            coordinates_out = zeros((n+1)*(m+1)*(p+1),3);
            connections_out = zeros(n*m*p,8);
            count = 1;
            for k = 1:p+1
                for j = 1:n+1
                    for i = 1:m+1
                        coordinates_out(count,:) = [(i-1)*x,(j-1)*y,(k-1)*z];
                        count = count + 1;
                    end
                end
            end
            count = 1;
            for k = 1:p
                for j = 1:n
                    for i = 1:m
                        connections_out(count,1) = 1+(i-1)+(j-1)*(m+1)+(k-1)*(m+1)*(n+1);
                        connections_out(count,2) = 1+i+(j-1)*(m+1)+(k-1)*(m+1)*(n+1);
                        connections_out(count,3) = i+(j)*(m+1)+(k-1)*(m+1)*(n+1);
                        connections_out(count,4) = 1+i+(j)*(m+1)+(k-1)*(m+1)*(n+1);
                        connections_out(count,5) = 1+(i-1)+(j-1)*(m+1)+k*(m+1)*(n+1);
                        connections_out(count,6) = 1+i+(j-1)*(m+1)+k*(m+1)*(n+1);
                        connections_out(count,7) = i+(j)*(m+1)+k*(m+1)*(n+1);
                        connections_out(count,8) = 1+i+(j)*(m+1)+k*(m+1)*(n+1);
                        count = count + 1;
                    end
                end
            end
        end
    end
end

    