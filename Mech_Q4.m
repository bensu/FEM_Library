classdef Mech_Q4 < Q4 & Mechanics
    methods
        function ele_out = Mech_Q4(nodes_in,material_in)
            ele_out = ele_out@Q4(nodes_in,material_in);
        end
    end    
end