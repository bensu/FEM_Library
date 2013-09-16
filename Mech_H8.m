classdef Mech_H8 < H8 & Mechanics
    methods
        function ele_out = Mech_H8(nodes_in,material_in)
            n_element_dofs = 0;
            n_node_dofs = 3;
            ele_out = ele_out@H8(n_node_dofs,n_element_dofs, ...
                                nodes_in,material_in);
        end
    end    
end