classdef ShellQ4 < Q4
    properties (Dependent)
        v3
        ele_thickness
    end
    methods
        function obj = ShellQ4(id_in,nodes,material,dofs_per_node,dofs_per_ele)
            if length(nodes) == 8
                nodes_in = Node.empty(length(nodes)/2,0);
                for n = 1:length(nodes)/2
                    xm = (nodes(n).get('coordinates') + nodes(n+4).get('coordinates'))/2;
                    v3 = (nodes(n).get('coordinates') - nodes(n+4).get('coordinates'))/2;
                    nodes_in(n) = Node(nodes(n).get('id'),xm);
                    nodes_in(n).set('v3',v3);
                end
            elseif length(nodes) == 4
                nodes_in = nodes;
            else error('Wrong Length of Nodes')
            end
            require(~isempty(nodes_in(1).get('v3')),'Nodes should have orientations')
            obj = obj@Q4(id_in,dofs_per_node,dofs_per_ele,nodes_in,material);
        end
        
        
        %% Plate-specific functions
        function u_out = u_matrix(element)
            xm = element.coordinates;
            V1 = zeros(3,element.nnodes);
            node = 1;
            for i = [-1,1]
                for j = [-1,1]
                    V1(:,node) = element.dN_deta([i,j,0])*xm;
                    node = node + 1;
                end
            end
            u_out = zeros(3,2,element.nnodes);
            for node = 1:element.nnodes
                u_out(:,2,node) = V1(:,node)/norm(V1(:,node)); 
                aux = cross(element.v3(:,node),u_out(:,2,node));
                u_out(:,1,node) = -aux/norm(aux);
            end
        end
        function t_at_node_out = t_at_node(element)
            t_at_node_out = zeros(element.nnodes,1);
            for node = 1:element.nnodes
                t_at_node_out(node) = norm(element.v3(:,node));
            end
        end  
        %% Setters & Getters
        function V3_out = get.v3(element)
            total_nodes = element.nnodes;
            V3_out = zeros(3,total_nodes);
            node_list = element.get('nodes');
            for i = 1:total_nodes
                V3_out(:,i) = node_list(i).get('v3');
            end
        end
        function t = get.ele_thickness(element)
            % t = get.ele_thickness(element)
            % Gives an approximation of the element thickness withs the
            % average of all the node thicknesses
            t = mean(element.t_at_node);
        end
            
        %% Geometry Functions
        function N_out = N(element,local_coords)
            [xi, eta, mu] = element.components(local_coords);
            require(isscalar(mu),'xi, eta, and mu should be scalar')
            require(-1<=mu && mu<=1,'mu should is not -1<=mu<=1')
            N_out = element.N@Q4([xi,eta]);
        end
        function dN_out = dN_dxi(element,local_coords)
            [xi, eta, mu] = element.components(local_coords);
            require(isscalar(mu),'xi, eta, and mu should be scalar')
            require(-1<=mu && mu<=1,strcat('mu should is not -1<=mu<=1: ',num2str(mu)))
            dN_out = element.dN_dxi@Q4([xi,eta]);
        end
        function dN_out = dN_deta(element,local_coords)
            [xi, eta, mu] = element.components(local_coords);
            require(isscalar(mu),'xi, eta, and mu should be scalar')
            require(-1<=mu && mu<=1,'mu should is not -1<=mu<=1')
            dN_out = element.dN_deta@Q4([xi,eta]);
        end
        function dN_out = dN_dmu(element,local_coords)
            [~, ~, mu] = element.components(local_coords);
            require(isscalar(mu),'xi, eta, and mu should be scalar')
            require(-1<=mu && mu<=1,'mu should is not -1<=mu<=1')
            dN_out = element.N(local_coords)*mu*element.v3;
        end
        function jacobian_out = jacobian(element,local_coords)
            [~, ~, mu] = element.components(local_coords);
            dNdxi = element.dN_dxi(local_coords);
            dNdeta = element.dN_deta(local_coords);
            N = element.N(local_coords);
            jacobian_out = [dNdxi;dNdeta;zeros(1,4)]*element.coordinates ...
                + [mu*dNdxi;mu*dNdeta;N]*element.v3'/2; % Cook 12.5-4 pag 360
            require(~near(det(jacobian_out),0),'detJac==0');
        end
    end
    
end