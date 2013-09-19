classdef Compound_Function < hgsetget
    properties
        element_component
        node_component
    end
    properties (Dependent)
        node_function
        element_function
    end
    methods
        function obj = Compound_Function(filler,total_number_of_nodes_in,dofs_per_node_in, ...
                            total_number_of_elements_in, dofs_per_element_in)
            node_component_in = Vector_Function(filler,total_number_of_nodes_in,dofs_per_node_in);
            element_component_in =  Vector_Function(filler,total_number_of_elements_in,dofs_per_element_in);
            obj.set('element_component',element_component_in);
            obj.set('node_component',node_component_in);
        end
        
        %% I/O
        function dofs = all_dofs(compound)
            dofs = [compound.get('node_component').all_dofs;
                       compound.get('element_component').all_dofs];
            
        end
        function dof_list_in(compound,dofs)
            nodes_dofs = dofs(1:compound.number_of_nodes_dofs);
            element_dofs = dofs((compound.number_of_nodes_dofs+1):end);
            assert(length(element_dofs)==compound.number_of_elements_dofs);
            compound.get('node_component').dof_list_in(nodes_dofs);
            compound.get('element_component').dof_list_in(element_dofs);
        end
        
        function fun_out = get.node_function(compound)
            fun_out = compound.get('node_component').get('component_function');
        end
        function fun_out = get.element_function(compound)
            fun_out = compound.get('element_component').get('component_function');
        end
        
        function add_to_node_function(compound,function_in)
        	compound.get('node_component').add_function(function_in);
        end
        function add_to_element_function(compound,function_in)
        	compound.get('element_component').add_function(function_in);
        end
        
        
        function clear(compound)
            compound.get('node_component').clear
            compound.get('element_component').clear
        end
        
        %% 
        
        function val_out = max_dof_value(compound)
            val_out = max(compound.all_dofs);
        end
        
        %% Helpers
        function num = number_of_elements(compound)
            num = compound.get('element_component').get('total_number_of_components');            
        end
        function num = number_of_nodes(compound)
            num = compound.get('node_component').get('total_number_of_components');            
        end
        function num = dofs_per_node(compound)
            num = compound.get('node_component').get('dofs_per_component');            
        end
        function num = dofs_per_element(compound)
            num = compound.get('node_element').get('dofs_per_component');           
        end
        function num = number_of_nodes_dofs(compound)
            num = compound.get('node_component').number_of_dofs;            
        end
        function num = number_of_elements_dofs(compound)
            num = compound.get('element_component').number_of_dofs;            
        end
    end                    
end