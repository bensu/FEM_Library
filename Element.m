classdef Element < hgsetget
    properties 
        nodes
        material
        parent
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

        function [Jac,Ndev] = jacob(obj,xi,eta,mu) %computes a 3x3 jacobian for an element in [xi,eta,mu] where X is it's coordinates
            checkv = [xi,eta,mu];
            for i = checkv
                if or(i>1,i<-1)
                    error('wrong args for xi')
                end
            end
            Ndev = Element.shfndev(xi,eta,mu);
            Jac = Ndev*obj.coordinates();
            if det(Jac)<0
                %error('det(Jacobian) < 0')
            end
        end
        
        function B = Bmatrix(obj,xi,eta,mu)
            AUX = zeros(6,9);
            AUX([1 10 18 22 26 35 42 47 51]) = 1;
            Tinv = inv(obj.jacob(xi,eta,mu));
            Ndevsparse = Element.shfndevsparse(xi,eta,mu);
            AUX2 = zeros(9);
            for i = 1:3
                AUX2((1+(i-1)*3:3*i),(1+(i-1)*3:3*i)) = Tinv;
            end
            B = AUX*AUX2*Ndevsparse;
        end
        
        function kele = kel(obj,gaussn)     %element stiffness matrix with gaussn order rule
            elmaterial = obj.get('material');
            kele = zeros(24);
            [gaussp gaussw] = lgwt(gaussn,-1,1);
            for i = 1:gaussn
                xi = gaussp(i);
                for j = 1:gaussn
                    eta = gaussp(j);
                    for k = 1:gaussn
                        mu = gaussp(k);
                        weight = gaussw(i)*gaussw(j)*gaussw(k);
                        B = obj.Bmatrix(xi,eta,mu);
                        kele = kele + weight*B'*elmaterial.get('C')*B*det(obj.jacob(xi,eta,mu));
                    end
                end
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

        
    end
    methods (Static)
        function N = shfn(xi,eta,mu)
            %returns a vector with shape functions
            N = zeros(1,8);
            count = 1;
            for k = [-1 1]
                for j = [-1 1]
                    for i = [-1 1]
                        N(count) = (1 + i*xi)*(1 + j*eta)*(1 + k*mu);
                        count = count + 1;
                    end
                end
            end
            
            N = N/8;
            
        end
        function Ndev = shfndev(xi,eta,mu)
            Ndev = [];
            for k = [-1 1]
                for j = [-1 1]
                    for i = [-1 1]
                        aux0 = [i*(1+j*eta)*(1+k*mu) j*(1+i*xi)*(1+k*mu) k*(1+j*eta)*(1+i*xi)]';
                        Ndev = [Ndev aux0];
                    end
                end
            end
            Ndev = Ndev/8;
        end
        function Ndevsparse = shfndevsparse(xi,eta,mu)
            AUX = Element.shfndev(xi,eta,mu);
            Ndevsparse = [];
            for i = 1:8
                aux0 = AUX(:,i);
                aux1 = [aux0 zeros(3,2)];
                aux2 = [zeros(3,1) aux0 zeros(3,1)];
                aux3 = [zeros(3,2) aux0];
                aux0 = [aux1;aux2;aux3];
                Ndevsparse = [Ndevsparse aux0];
            end
        end
        function connected_nodes = connectednodes(nodenum)
            AUX = [2 3 5; 1 4 6; 1 4 7; 2 3 8; 1 6 7; 2 5 8; 3 5 8; 4 6 7];
            connected_nodes = AUX(nodenum,:);
        end
    end

end
