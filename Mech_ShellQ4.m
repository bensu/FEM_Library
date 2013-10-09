classdef Mech_ShellQ4 < Mech_Shell & ShellQ4
    properties (Constant)
        mech_node_dofs_id = 1:5;    % Dofs Ids Mechanics takes for each node
        mech_ele_dofs_id = [];       % Dofs Ids Mechanics takes for each element
        n_stress_components = 6;
    end
    methods
        function obj = Mech_ShellQ4(id_in,nodes,material)
            require(length(nodes)==4 || length(nodes) == 8,'Wrong node_list size');
            obj = obj@ShellQ4(id_in,nodes,material,5,0);
        end
    end
end