classdef Mech_H8 < H8 & Mechanics
    properties (Constant)
        mech_node_dofs_id = 1:3;    % Dofs Ids Mechanics takes for each node
        mech_ele_dofs_id = [];       % Dofs Ids Mechanics takes for each element
        n_stress_components = 6;
    end
    methods
        function ele_out = Mech_H8(id_in,nodes_in,material_in)
            n_element_dofs = 0;
            n_node_dofs = 3;
            ele_out = ele_out@H8(id_in,n_node_dofs,n_element_dofs, ...
                                nodes_in,material_in);
        end
    end    
end