classdef Reactions < Compound_Function
    methods 
        function react_out = Reactions(total_number_of_nodes_in,dofs_per_node_in, ...
                            total_number_of_elements_in, dofs_per_element_in)
            react_out = react_out@Compound_Function(0,total_number_of_nodes_in,dofs_per_node_in, ...
                            total_number_of_elements_in, dofs_per_element_in);    
        end
    end
end