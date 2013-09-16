classdef (Abstract) Element < hgsetget
    properties
        nodes
        material
        parent
        n_ele_dofs
        n_dofs_per_node
    end
    methods (Abstract)
        N_out = N(element,xi,eta,mu)
        dN_out = dN_dxi(element,xi,eta,mu)
        dN_out = dN_deta(element,xi,eta,mu)
        dN_out = dN_dmu(element,xi,eta,mu)
        dN_out = DN(element,xi,eta,mu)
    end
        
    methods
        function obj = Element(nodes,material)
            set(obj,'nodes',nodes);
            set(obj,'material',material);
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
        function coordinat = coordinates(obj)
            nodelist = obj.get('nodes');
            coordinat = zeros(numel(nodelist),3);
            for i = 1:numel(nodelist)
                nodei = nodelist(i);
                coordinat(i,:) = nodei.get('coordinates');
            end
        end

        function jac_out = jacobian(element,xi,eta,mu) 
            % jac_out [3x3] = jacobian(obj,xi,eta,mu) 
            checkv = [xi,eta,mu];
            for i = checkv
                if or(i>1,i<-1)
                    error('wrong args for xi,eta, or mu')
                end
            end
            jac_out = element.DN(xi,eta,mu)*element.coordinates();
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
        function n_out = n_nodes(element)
            % n_out = n_nodes(element)
            % Number of nodes in the element
            n_out = length(element.get('nodes'));
        end

        
    end
 

end
