classdef ME_ShellQ4 <  Mech_Shell & ShellQ4
    properties (Constant)
        mech_node_dofs_id = 1:5;     % Dofs Ids Mechanics takes for each node
        mech_ele_dofs_id = [];       % Dofs Ids Mechanics takes for each element
        electric_node_dofs_id = [];
        electric_ele_dofs_id = 1:2;
        n_stress_components = 7;     % Six tensions [N/m2] and a Dielectric Displacement Dzz
    end
    methods
        function obj = ME_ShellQ4(id_in,nodes,material)
            require(length(nodes)==4 || length(nodes) == 8,'Wrong node_list size');
            obj = obj@ShellQ4(id_in,nodes,material,5,2);
        end
        function Ndevs_out = DNsparse(element,local_coords)
            Ndevs_out = blkdiag(element.DNsparse@Mech_Shell(local_coords), ...
                                                                [-1 1]);
        end 
        function H_out = dU_to_strain(element)
            H_out = blkdiag(element.dU_to_strain@Mechanics,-1);
        end
        function invT_out = invT(element,local_coords)
            invT_out = blkdiag(element.invT@Mechanics(local_coords), ...
                                    1/element.ele_thickness);
        end
        function C_out = C(element)
            C_out = blkdiag(element.C@Mech_Shell,-element.material.get('Permitivity3'));
            C_out([7 14 43 44]) = -element.material.get('Piezo13');
        end
    end
end