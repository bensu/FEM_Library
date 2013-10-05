classdef Mech_Q4 < Q4 & Mechanics
    properties (Constant)
        mech_node_dofs_id = 1:2;    % Dofs Ids Mechanics takes for each node
        mech_ele_dofs_id = [];       % Dofs Ids Mechanics takes for each element
        n_stress_components = 3;    % Indicates Plain Strain
    end
    methods
        function ele_out = Mech_Q4(id_in,nodes_in,material_in)
            n_element_dofs = 0;
            n_node_dofs = 2;        % u, v
            ele_out = ele_out@Q4(id_in,n_node_dofs,n_element_dofs,...
                            nodes_in,material_in);
        end
    end    
end