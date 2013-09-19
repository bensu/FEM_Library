classdef BC < Compound_Function
    methods
        function bc_out = BC(total_number_of_nodes_in,dofs_per_node_in, ...
                            total_number_of_elements_in, dofs_per_element_in)
            bc_out = bc_out@Compound_Function(true,total_number_of_nodes_in,dofs_per_node_in, ...
                            total_number_of_elements_in, dofs_per_element_in);    
        end
        function fix_node_by_id(bc,node_id)
            bc.get('node_component').edit_component_by_id(node_id, ...
            	false(1,bc.get('node_component').get('dofs_per_component')));
        end
        function support_node_by_id(bc,node_id,direction)
            bc.get('node_component').edit_component_part_by_id(...
                                    node_id,direction,false);
        end
        function support_element_by_id(bc,ele_id,direction)
            bc.get('element_component').edit_component_part_by_id(...
                                    ele_id,direction,false);
        end
        function support_all_elements_in_dir(bc,direction)
            for ele_id = 1:bc.number_of_elements
                bc.support_element_by_id(ele_id,direction)
            end
        end
        function obj = fixed_face(obj,new_face)
            previous = obj.node_function();
            nodelist = new_face.get('nodelist');
            for i = 1:length(nodelist)
                previous(nodelist(i),:) = false(1,3);
            end
            obj.set('nodelist',VectorField.xyzin(previous));
        end
        function face_constraint(obj,face,coordnum,mesh)
            nodelistids = mesh.facenodelist(face);
            for i = 1:length(nodelistids)
                obj.support_node_by_id(nodelistids(i),coordnum);
            end
        end
        function simply_supported_face(bc,face,coordnum,mesh) %missing a point
            coords = mesh.get('coordinates');
            facelist = mesh.facenodelist(face);
            %first fixed point
            bc.fix_node_by_id(facelist(1));
            %face
            bc.face_constraint(face,coordnum,mesh);
            
            dim = size(coords,2);
            switch size(coords,2)
                case 2
%                     orthogonal_coord = mod(coordnum,2)+1;
                case 3
                    v1 = zeros(1,dim);
                    v1(coordnum) = 1;
                    v2 = coords(facelist(1),:) - coords(facelist(2),:);
                    orthogonal_coord = find(cross(v1,v2));
                    bc.support_node_by_id(facelist(2),orthogonal_coord); 
            end
                                
                     
        end  
        function plot(obj,coordinates,color)
            hold on
            xyz = ~obj.node_function();
            quiver3(coordinates(:,1),coordinates(:,2),coordinates(:,3),xyz(:,1),xyz(:,2),xyz(:,3),color)
        end
    end
end
    