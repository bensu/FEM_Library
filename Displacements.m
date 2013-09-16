classdef Displacements < Compound_Function
    methods 
        function dis_out = Displacements(total_number_of_nodes_in,dofs_per_node_in, ...
                            total_number_of_elements_in, dofs_per_element_in)
            dis_out = dis_out@Compound_Function(0,total_number_of_nodes_in,dofs_per_node_in, ...
                            total_number_of_elements_in, dofs_per_element_in);    
        end
    end
end