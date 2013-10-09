classdef (Abstract) Mechanics < Physics
    properties (Abstract,Constant)
        mech_node_dofs_id    % Dofs Ids Mechanics takes for each node
        mech_ele_dofs_id     % Dofs Ids Mechanics takes for each element
        n_stress_components
    end
    properties (Dependent)
        dis
        mech_dis
        mech_n_node_dofs
        mech_n_ele_dofs
        mech_n_dofs
        mech_dofs_id
        n_dofs_per_node %?????
    end
    methods
        function n_out = get.n_dofs_per_node(element)
            n_out = element.dim;
        end
        function Tinv = invT(element,local_coords)
            dim_aux = element.dim;
            Tinv = zeros(dim_aux^2);
            invJ = inv(element.jacobian(local_coords));
            for i = 1:dim_aux
                range = (1+(i-1)*dim_aux:dim_aux*i);
                Tinv(range,range) = invJ;
            end
        end
        function stress_voight = node_stress(element)
            % stress_voigth = node_stress(element)
            % dis = [n_nodes*3 x 1]
            % stress_voigth [n_nodes xn_stress_components]
            % Takes the nodal displacements and returns the stress tensor in
            % Voight vector notation: [sigma_xx s_yy s_zz tau_xy t_xz t_yz];
            % returns one column per node.
            stress_voight = (element.C*element.node_strain')';
        end
        function strain_voight = node_strain(element)
            % strain_voigth = stress_from_dis(element,dis)
            % dis = [n_nodes*3 x 1]
            % strain_voigth [element.nnodes x n_stress_components]
            % Takes the nodal displacements and returns the stress tensor in
            % Voight vector notation: [sigma_xx s_yy s_zz tau_xy t_xz t_yz]';
            strain_voight = zeros(element.nnodes,element.n_stress_components);
            for n = 1:element.nnodes
                l_coords = element.node_local_coords(n,:);
                strain_voight(n,:) = (element.B(l_coords)*element.dis)';
            end
        end
        function B_out = B(element,local_coords)
            % B_out [stress_components x element_dofs] = B(element,xi,eta,mu)
            % Computes the B matrix generically for any element in
            % structural mechanics
            
            % NOTE - Needs refactoring, since now its specifically
            % programmed for H8

            B_out = element.dU_to_strain*element.invT(local_coords)* ...
                        element.DNsparse(local_coords);
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
            % programmed for H8 & Q4
            K_out = zeros(element.n_dofs);
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
        function H_out = dU_to_strain(element)
            AUX = zeros(6,9);
            AUX([1 10 18 22 26 35 42 47 51]) = 1;
            if element.dim == 3
                H_out = AUX;
            elseif element.dim == 2
                H_out = AUX([1 2 4],[1 2 4 5]);
            end
        end
        
        %% Setters & Getters
        
        function dis_out = get.mech_dis(element)
            dis_out = element.dis(element.mech_dofs_id);
        end
        function dof_list = get.mech_dofs_id(element)
            ind = [];
            for n = 1:element.nnodes
                ind = [ind; index_range(element.dofs_per_node,n,element.mech_node_dofs_id)];
            end
            dof_list = [ind; element.n_node_dofs + element.mech_ele_dofs_id];
        end
       
        function n_dof = get.mech_n_node_dofs(element)
            n_dof = length(element.mech_node_dofs_id);
        end
        function n_dof = get.mech_n_ele_dofs(element)
            n_dof = length(element.mech_ele_dofs_id);
        end
        function n_dof = get.mech_n_dofs(element)
            n_dof = element.mech_n_node_dofs + element.mech_n_ele_dofs;
        end
        function dis = get.dis(element)
            dis = element.get('dof_dis');
        end
    end
            
        
end

