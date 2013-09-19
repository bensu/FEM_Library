classdef Mesh < hgsetget
    properties
        element_type
        coordinates
        connections
        material
        outer_faces
        nodesout
        parent
        ncambio
        n_node_dofs         % Dofs per node
        n_element_dofs      % Dofs that belong only to the element, not the nodes
    end
    
    methods (Static)
        function L = NXorder()
            L = [6 2 5 1 7 3 8 4];
        end
        function sca = tometer()
            sca = 1e3;
        end
        function obj = mesh_import(coordstr,connectstr,mate)
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
            obj = Mesh(coordinates_in,connections_in,mate);
        end
        function obj = mesh_import2(coordstr,connectstr,mate)
            connections_in = load(connectstr);
            coordinates_in = load(coordstr);
            connections2 = zeros(size(connections_in));
            for i = 1:size(coordinates_in,1)
                nodenum = coordinates_in(i,1);
                connections2(connections_in==nodenum) = i;
            end
            
            coordinates_in = coordinates_in/1e4;
            coordinates_in = coordinates_in/Mesh.tometer();
            obj = Mesh(coordinates_in,connections_in,mate);
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
    
    methods

        function obj = Mesh(type,coordinates,connections,material)
            set(obj,'coordinates',coordinates);
            set(obj,'connections',connections);
            set(obj,'material',material);
            set(obj,'element_type',type);
            ele = obj.element_create(1);
            set(obj,'n_node_dofs',ele.get('n_node_dofs'));       % Dofs per node
            set(obj,'n_element_dofs',ele.get('n_element_dofs'));  
        end

        %% Coordinates Methods
        
        function value = mincoordval(obj,coordnum)
            coordinatesl = obj.get('coordinates');
            [~, ind] = min(abs(coordinatesl(:,coordnum)));
            value = coordinatesl(ind,coordnum);
        end
        function value = maxcoordval(obj,coordnum)
            coordinatesl = obj.get('coordinates');
            [~, ind] = max(abs(coordinatesl(:,coordnum)));
            value = coordinatesl(ind,coordnum);
        end
        function shift(obj,dis_vector)
            coords = obj.get('coordinates');
            coords2 = zeros(size(coords));
            for n = 1:obj.nnodes()
                coords2(n,:) = coords(n,:) + dis_vector;
            end
            obj.set('coordinates',coords2)
        end
        
        %% Face Methods
        
        function add_outer_face(obj,newface)
            obj.set('outer_faces',[obj.get('outer_faces');newface])
        end
        function outfaces = outface(mesh)
            % outfaces = outface(mesh)
            % Finds the outer nodes and returns an out face object
            connect = mesh.get('connections');
            nodelist = [];
            par_plot = [];
            for nodenum = 1:mesh.nnodes()   % Loop through all the nodes
                elelist = mesh.elementsofnode(nodenum);
                % Finds the elements the node belongs too
                element = mesh.element_create(elelist(1));
                if length(elelist) < element.n_nodes    % If
                    nodelist = [nodelist nodenum];
                    for i = 1:length(elelist)
                    	allnodes = connect(elelist(i),:);  %should be [n_nodesx1]
                        other_nodes = allnodes(element.connected_nodes(nodenum == allnodes));
                        for j = 1:length(other_nodes)
                            if any(other_nodes(j)==nodelist)
                                par_plot = [par_plot; nodenum other_nodes(j)];
                            end
                        end
                    end
                end
            end
%             big_elelist = unique(big_elelist);
            outfaces = Face(nodelist,mesh);
            outfaces.set('par_plot',par_plot);
        end
        function facenodelist = facenodelist(mesh,face)
            facenodelist = face.get('nodelist');
        end 
        %% Node Methods
        function innodes = innernodes(obj)
            outface = obj.outface();
            outnodelist = obj.facenodelist(outface);
            innodes = [];
            count = 1;
            for i = 1:obj.nnodes()
                if ~any(outnodelist==i)
                    innodes(count) = i;
                    count = count + 1;
                end
            end
        end
        function nodesin = nodes_in(obj)
            nodesin = 1:size(obj.get('coordinates'),1);      
            nodesin(obj.get('nodesout')) = [];   
            ncambio_aux = obj.get('ncambio');
            nodesout = reshape(ncambio_aux(nodesin),1,[]); % WHAT IS HERE?
        end
        function nodesout = nodes_out(obj)
            ncambio_aux = obj.get('ncambio');
            nodesout = reshape(ncambio_aux(obj.get('nodesout')),1,[]);      
        end 
        function ind = nodesindex(obj,nodelist)
            ndofpernode = obj.ndofspernode();
            ind = zeros(length(nodelist),1);
            for i = 1:length(nodelist)
                range = (ndofpernode*(i-1)+1):(ndofpernode*i);
                ind(range) = (ndofpernode*(nodelist(i)-1)+1):(ndofpernode*nodelist(i));
            end
        end
        function elementlist = elementsofnode(obj,nodenum)
            connections1  = obj.get('connections');
            elementlist = [];
            for i = 1:size(connections1,1)
                if any(nodenum==connections1(i,:))
                    elementlist = [elementlist i];
                end
            end
        end
        function randominner(obj,delta)
            innodes = obj.innernodes();
            if ~isempty(innodes)
                coord = obj.get('coordinates');
                coord(innodes(1),1) = coord(innodes(1),1) + delta;
                obj.set('coordinates',coord);
            end
        end
        function nodenum = findnode(obj,x0)
            nodenum = [];
            coordinates_aux = obj.get('coordinates');
            tol = 5e-3;
            for i = 1:obj.nnodes()
                if norm(coordinates_aux(i,:)-x0)<tol;
                    nodenum = i;
                end
            end
        end
        
        function elementcoord_out = elementcoord(obj,ele)
            new_ele = obj.elementcreate(ele);
            elementcoord_out = median(new_ele.coordinates());
        end
        
        %% Assembly Methods
        function element = element_create(mesh,ele)
            connect = mesh.get('connections');
            coordinat = mesh.get('coordinates');
            for i = 1:mesh.nodesperelement()
                nodos(i) = Node(connect(ele,i),coordinat(connect(ele,i),:),mesh);
            end
            % BIG PROBLEM HERE -> How to initialize by element type?
            switch mesh.element_type
                case 'Mech_H8'
                    element = Mech_H8(nodos,mesh.get('material'));   % HARDCODED
                case 'Mech_Q4'
                    element = Mech_Q4(nodos,mesh.get('material'));
                otherwise
                    error(strcat('Element Type: ',mesh.element_type,' not found'));
            end
        end
        function stiffness = stiff(obj)
                nnel = obj.nnel();
                stiffness = sparse(obj.sdof(),obj.sdof());
                gaussn = 2;
                for ele = 1:nnel
                    element = obj.element_create(ele);
                    kele = element.K(gaussn);
                    ind = element.ldofsid();
                    stiffness(ind,ind) = stiffness(ind,ind) + kele;
                end
        end


        %% Property Methods
        function out = nnodes(obj)          %total number of nodes
            out = size(get(obj,'coordinates'),1);
        end
        function out = nnel(obj)            %total number of elements
            out = size(get(obj,'connections'),1);
        end
        function out = ndofspernode(obj)    %number of dofs per node
            out = size(get(obj,'coordinates'),2);
        end
        function out = sdof(obj)            %total system dofs
            out = obj.nnodes()*obj.ndofspernode();
        end
        function out = nodesperelement(obj) %nodes per element
            out = size(get(obj,'connections'),2);
        end
        function out = dim(obj)    %number of dofs per node
            out = size(get(obj,'coordinates'),2);
        end
        
        
        %% Plot
        
        function plot(mesh,colour)
            hold on
            npoints = 2;
            outfaces = mesh.outface();
            coordinates_aux = mesh.get('coordinates');
            par_plot_now = outfaces.get('par_plot');
            X = zeros(npoints,3);
            for i = 1:size(par_plot_now,1)
                for j = 1:mesh.dim
                    X(:,j) = linspace(coordinates_aux(par_plot_now(i,1),j), ...
                               coordinates_aux(par_plot_now(i,2),j),npoints)';
                end
                plot3(X(:,1),X(:,2),X(:,3),colour);
            end
            hold off
        end
        
                %% Plot
        
        function plot_function(mesh,vf,dim,type,color)
        % plot(vf,dim,type,color)
        % dim: Which dimensions of the Vector functio to Plot
        %      If it is a displacement, U = [u,v,w], u can be plotted by
        %      specifying dim = 1.
        % type: patch for 1 dimension
        %       quiver for 2 or 3 dimensions 
            coords = mesh.get('coordinates');
            if mesh.dim == 2
                coords = [coords zeros(size(coords,1),1)];
            end
            switch type
                case 'quiver'
                    require(any(vf.dofs_per_node == [2,3]), ...
                                    'Wrong Dimension for Quiver');
                    quiver3(coords(:,1),coords(:,2),coords(:,3), ...
                        vf.node_function(:,1),vf.node_function(:,2), ...
                            vf.node_function(:,3));
                case 'patch'
                    require(length(dim)==1 && vf.dofs_per_component==1)
                    
            end        
        end
    end
end

    