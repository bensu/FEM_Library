classdef (Abstract) H8 < Element
    methods
        function ele_out = H8(nodes_in,material_in)
            ele_out = ele_out@Element(nodes_in,material_in);
        end
        function N_out = N(element,xi,eta,mu)
        % N_out [1x8] = N(element,xi,eta,mu) 
        % Linear H8 Shape Functions
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
        function dN_out = dN_dxi(element,xi,eta,mu)
        % dN_out [1x8] = dN_dxi(element,xi,eta,mu)
        % Linear H8 Shape Functions differenciated by xi. Do not depend on xi.
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
        function dN_out = dN_deta(element,xi,eta,mu)
        % dN_out [1x8] = dN_deta(element,xi,eta,mu)
        % Linear H8 Shape Functions differenciated by eta. Do not depend on eta.
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
        function dN_out = dN_dmu(element,xi,eta,mu)
        % dN_out [1x8] = dN_dmu(element,xi,eta,mu)
        % Linear H8 Shape Functions differenciated by mu. Do not depend on mu.
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
        
        
        function Ndev = DN(element,xi,eta,mu)
            % Ndev [3x8] = DN(element,xi,eta,mu)
            % Groups all the derivatives of the shape functions in a matrix
            Ndev = [element.dN_dxi(xi,eta,mu); 
                   element.dN_deta(xi,eta,mu);
                   element.dN_dmu(xi,eta,mu)];
        end
        function connected_nodes = connected_nodes(element,nodenum)
            AUX = [2 3 5; 1 4 6; 1 4 7; 2 3 8; 1 6 7; 2 5 8; 3 5 8; 4 6 7];
            connected_nodes = AUX(nodenum,:);
        end
    end
end
