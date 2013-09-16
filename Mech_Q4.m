classdef Mech_Q4 < Q4 & Mechanics
    methods
        function ele_out = Mech_Q4(nodes_in,material_in)
            n_element_dofs = 0;
            n_node_dofs = 2;        % u, v
            ele_out = ele_out@Q4(n_node_dofs,n_element_dofs,...
                            nodes_in,material_in);
        end
    end    
end