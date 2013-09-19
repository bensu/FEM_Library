classdef (Abstract) Mechanics < hgsetget
    properties (Dependent)
        n_stress_components
        dU_to_strain
    end
    methods
        function Tinv = invT(element,local_coords)
            dim_aux = element.dim;
            Tinv = zeros(dim_aux^2);
            invJ = inv(element.jacobian(local_coords));
            for i = 1:dim_aux
                range = (1+(i-1)*dim_aux:dim_aux*i);
                Tinv(range,range) = invJ;
            end
        end
        function B_out = B(element,local_coords)
            % B_out [stress_components x element_dofs] = B(element,xi,eta,mu)
            % Computes the B matrix generically for any element in
            % structural mechanics
            
            % NOTE - Needs refactoring, since now its specifically
            % programmed for H8

            Tinv = element.invT(local_coords);
            Ndevsparse = element.DNsparse(local_coords);
            B_out = element.dU_to_strain*Tinv*Ndevsparse;
        end   
        function C_out = C(element)
            % C_out [6x6] = C(element)
            % Returns the C matrix for an Isotropic Material obeying Hook's
            % Law
            E = element.material.get('Young_Modulus');
            nu = element.material.get('Poisson_Coefficient');
            C_3D = E/((1+nu)*(1-2*nu))*[1-nu nu   nu   0          0          0          
                       nu   1-nu nu   0          0          0  
                       nu   nu   1-nu 0          0          0
                       0    0    0    (1-2*nu)/2 0          0
                       0    0    0    0          (1-2*nu)/2 0
                       0    0    0    0          0          (1-2*nu)/2];
            if element.dim == 3
                C_out = C_3D;
            elseif element.dim == 2
                C_out = C_3D([1 2 4],[1 2 4]);
            end
        end
        function K_out = K(element,gauss_order)
            % K_out [ele_dofs x element_dofs] = K_dofs(element,gauss_order,xi,eta,mu)
            % Computes the stiffness matrix generically for any element in
            % structural mechanics
            
            % NOTE - Needs refactoring, since now its specifically
            % programmed for H8
            K_out = zeros(element.dim*element.n_nodes);
            [gaussp, gaussw] = lgwt(gauss_order,-1,1);
            for i = 1:gauss_order
                xi = gaussp(i);
                for j = 1:gauss_order
                    eta = gaussp(j);
                    if element.dim == 3
                        for k = 1:gauss_order
                            mu = gaussp(k);
                            weight = gaussw(i)*gaussw(j)*gaussw(k);
                            B = element.B([xi,eta,mu]);
                            K_out = K_out + weight*B'*element.C*B* ...
                                        det(element.jacobian([xi,eta,mu]));
                        end
                    elseif element.dim == 2
                            weight = gaussw(i)*gaussw(j);
                            B = element.B([xi,eta]);
                            K_out = K_out + weight*B'*element.C*B* ...
                                        det(element.jacobian([xi,eta]));
                    end
                        
                end
            end
        end
        function Ndevsparse = DNsparse(element,local_coords)
            % Ndevsparse = DNsparse(element,xi,eta,mu)
            % Creates an Auxiliary Matrix for B
            AUX = element.DN(local_coords);
            Ndevsparse = [];
            for i = 1:element.n_nodes
                if element.dim == 2
                    aux0 = blkdiag(AUX(:,i),AUX(:,i));
                elseif element.dim == 3
                    aux0 = blkdiag(AUX(:,i),AUX(:,i),AUX(:,i));
                end
                Ndevsparse = [Ndevsparse aux0];
            end
        end
        
        %% Setters & Getters
        function n_sigma = get.n_stress_components(element)
            if element.dim == 2
                n_sigma = 3;
            elseif element.dim == 3
                n_sigma = 6;
            end
        end
        
        function H_out = get.dU_to_strain(element)
            AUX = zeros(6,9);
            AUX([1 10 18 22 26 35 42 47 51]) = 1;
            if element.dim == 3
                H_out = AUX;
            elseif element.dim == 2
                H_out = AUX([1 2 4],[1 2 4 5]);
            end
        end
    end
            
        
end

