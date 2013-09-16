classdef Mech_H8 < H8 & Mechanics
    methods
        function ele_out = Mech_H8(nodes_in,material_in)
            ele_out = ele_out@H8(nodes_in,material_in);
        end
    end    
end