classdef (Abstract) Element < hgsetget
    properties
        nodes
        material
        parent
        n_element_dofs
        n_node_dofs
    end
    properties (Dependent)
        dim
        n_nodes
    end
    methods (Abstract)
        N_out = N(element,local_coords)
        dN_out = dN_dxi(element,local_coords)
        dN_out = dN_deta(element,local_coords)
        dN_out = dN_dmu(element,local_coords)
        dN_out = DN(element,local_coords)
    end
        
    methods
        function obj = Element(n_node_dofs_in,n_element_dofs_in,nodes, ...
                                        material)
            set(obj,'nodes',nodes);
            set(obj,'material',material);
            set(obj,'n_node_dofs',n_node_dofs_in);
            set(obj,'n_element_dofs',n_element_dofs_in);
        end
        function nodeids = nodeidlist(obj)
        	nodelist = obj.get('nodes');
            nodeids = zeros(size(nodelist));
            for i = 1:numel(nodelist)
                node = nodelist(i);
                nodeids(i) = node.get('id');
            end
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
            nodesids = obj.nodeidlist();
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
        
        
        
        function ind = ldofsid(obj)
            ind = [];
            nodelist = obj.get('nodes');
            for i = 1:length(nodelist)
                node = nodelist(i);
                ind = [ind node.ldofsid()];
            end            
        end

        %% Setters and Getters
        
        function [xi, eta, mu] = components(element,local_coords)
            if element.dim == 3 && length(local_coords) == 3
                xi = local_coords(1);
                eta = local_coords(2);
                mu = local_coords(3);
            elseif element.dim == 2 && length(local_coords) == 2
                xi = local_coords(1);
                eta = local_coords(2);
                mu = NaN;
            else error('Dimension missmatch local_coords and element dim')
            end
                
        end
        function dim_out = get.dim(element)
            node_list = element.get('nodes');
            dim_out = length(node_list(1).get('coordinates'));
        end
        function dim_out = get.n_nodes(element)
            dim_out = length(element.get('nodes'));
        end
    end
 

end
