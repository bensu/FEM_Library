classdef (Abstract) Mechanics < hgsetget
    properties (Constant)
        n_stress_components = 6;
    end
    methods
        function B_out = B(element,xi,eta,mu)
            % B_out [stress_components x element_dofs] = B(element,xi,eta,mu)
            % Computes the B matrix generically for any element in
            % structural mechanics
            
            % NOTE - Needs refactoring, since now its specifically
            % programmed for H8
            
            AUX = zeros(Mechanics.n_stress_components,9);
            AUX([1 10 18 22 26 35 42 47 51]) = 1;
            Tinv = inv(element.jacobian(xi,eta,mu));
            Ndevsparse = element.DNsparse(xi,eta,mu);
            AUX2 = zeros(9);
            for i = 1:3
                AUX2((1+(i-1)*3:3*i),(1+(i-1)*3:3*i)) = Tinv;
            end
            B_out = AUX*AUX2*Ndevsparse;
        end   
        function C_out = C(element)
            % C_out [6x6] = C(element)
            % Returns the C matrix for an Isotropic Material obeying Hook's
            % Law
            E = element.material.get('Young_Modulus');
            nu = element.material.get('Poisson_Coefficient');
            C_out = E/((1+nu)*(1-2*nu))*[1-nu nu   nu   0          0          0          
                       nu   1-nu nu   0          0          0  
                       nu   nu   1-nu 0          0          0
                       0    0    0    (1-2*nu)/2 0          0
                       0    0    0    0          (1-2*nu)/2 0
                       0    0    0    0          0          (1-2*nu)/2];
        end
        function K_out = K(element,gauss_order)
            % K_out [ele_dofs x element_dofs] = K_dofs(element,gauss_order,xi,eta,mu)
            % Computes the stiffness matrix generically for any element in
            % structural mechanics
            
            % NOTE - Needs refactoring, since now its specifically
            % programmed for H8
            K_out = zeros(24);
            [gaussp, gaussw] = lgwt(gauss_order,-1,1);
            for i = 1:gauss_order
                xi = gaussp(i);
                for j = 1:gauss_order
                    eta = gaussp(j);
                    for k = 1:gauss_order
                        mu = gaussp(k);
                        weight = gaussw(i)*gaussw(j)*gaussw(k);
                        B = element.B(xi,eta,mu);
                        K_out = K_out + weight*B'*element.C*B* ...
                                    det(element.jacobian(xi,eta,mu));
                    end
                end
            end
        end
        function Ndevsparse = DNsparse(element,xi,eta,mu)
            % Ndevsparse = DNsparse(element,xi,eta,mu)
            % Creates an Auxiliary Matrix for B
            AUX = element.DN(xi,eta,mu);
            Ndevsparse = [];
            for i = 1:element.n_nodes
                aux0 = AUX(:,i);
                aux1 = [aux0 zeros(3,2)];
                aux2 = [zeros(3,1) aux0 zeros(3,1)];
                aux3 = [zeros(3,2) aux0];
                aux0 = [aux1;aux2;aux3];
                Ndevsparse = [Ndevsparse aux0];
            end
        end
    end
end
