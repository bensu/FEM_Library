classdef (Abstract) H8 < Element
    properties (Constant)
        nnodes = 8;
        node_connectivity = [1 2;1 3;1 5;2 4;2 6;3 4;3 7;4 8;5 6;5 7;6 8;7 8];
        face_connectivity = [1 2 4 3;1 2 6 5;1 3 7 5;8 6 5 7;8 4 3 7;8 4 2 6];
        node_local_coords = [-1 -1 -1;-1 1 -1;1 -1 -1;1 1 -1;-1 -1 1;-1 1 1;1 -1 1;1 1 1];
    end
    methods
        function ele_out = H8(id_in,n_node_dofs,n_element_dofs,nodes_in,material_in)
            require(length(nodes_in(1).get('coordinates'))==3, ...
                'Mesh should be 3D');
            ele_out = ele_out@Element(id_in,n_node_dofs,n_element_dofs, ...
                                nodes_in,material_in);
        end
        function N_out = N(element,local_coords)
        % N_out [1x8] = N(element,local_coords) 
        % Linear H8 Shape Functions
            [xi, eta, mu] = element.components(local_coords);
            N_out = zeros(1,8);
            count = 1;
            for k = [-1 1]
                for j = [-1 1]
                    for i = [-1 1]
                        N_out(count) = (1 + i*xi)*(1 + j*eta)*(1 + k*mu);
                        count = count + 1;
                    end
                end
            end
            N_out = N_out/8;
        end
        function dN_out = dN_dxi(element,local_coords)
        % dN_out [1x8] = dN_dxi(element,xi,eta,mu)
        % Linear H8 Shape Functions differenciated by xi. Do not depend on xi.
            [~, eta, mu] = element.components(local_coords);    
            dN_out = zeros(1,8);
            count = 1;
            for k = [-1 1]
                for j = [-1 1]
                    for i = [-1 1]
                        dN_out(count) = i*(1+j*eta)*(1+k*mu);
                        count = count + 1;
                    end
                end
            end
            dN_out = dN_out/8;    
        end
        function dN_out = dN_deta(element,local_coords)
        % dN_out [1x8] = dN_deta(element,xi,eta,mu)
        % Linear H8 Shape Functions differenciated by eta. Do not depend on eta.
            [xi, ~, mu] = element.components(local_coords);    
            dN_out = zeros(1,8);
            count = 1;
            for k = [-1 1]
                for j = [-1 1]
                    for i = [-1 1]
                        dN_out(count) = j*(1+i*xi)*(1+k*mu);
                        count = count + 1;
                    end
                end
            end
            dN_out = dN_out/8;    
        end
        function dN_out = dN_dmu(element,local_coords)
        % dN_out [1x8] = dN_dmu(element,xi,eta,mu)
        % Linear H8 Shape Functions differenciated by mu. Do not depend on mu.
            [xi, eta, ~] = element.components(local_coords);   
            dN_out = zeros(1,8);
            count = 1;
            for k = [-1 1]
                for j = [-1 1]
                    for i = [-1 1]
                        dN_out(count) = k*(1+j*eta)*(1+i*xi);
                        count = count + 1;
                    end
                end
            end
            dN_out = dN_out/8;    
        end
        
        
        function Ndev = DN(element,local_coords)
            % Ndev [3x8] = DN(element,xi,eta,mu)
            % Groups all the derivatives of the shape functions in a matrix
            Ndev = [element.dN_dxi(local_coords); 
                   element.dN_deta(local_coords);
                   element.dN_dmu(local_coords)];
        end
        function connected_nodes = connected_nodes(element,nodenum)
            AUX = [2 3 5; 1 4 6; 1 4 7; 2 3 8; 1 6 7; 2 5 8; 3 5 8; 4 6 7];
            connected_nodes = AUX(nodenum,:);
        end
    end
end
