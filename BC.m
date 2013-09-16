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
        function simply_supported_face(obj,face,coordnum,mesh) %missing a point
            coordinates = mesh.get('coordinates');
            facelist = mesh.facenodelist(face);
            %first fixed point
            nodeid = facelist(1);
            obj.fix_node_by_id(nodeid);
            %face
            obj.face_constraint(face,coordnum,mesh);
            u = zeros(3,1);
            u(coordnum) = 1;
            iteration = 1;
            log = 0;
            iterationlimit = 1;
            while and(log == 0,iterationlimit<100)
                iteration = iteration + 1;
                coordnum2 = find(cross(u,coordinates(nodeid,:) - coordinates(nodeid+1,:)));
                if ~isempty(coordnum2)
                    obj.support_node_by_id(nodeid+1,coordnum2);
                    log = 1;
                else nodeid = nodeid + 1;
                end
            end
          
        end  
        function plot(obj,coordinates,color)
            hold on
            xyz = ~obj.node_function();
            quiver3(coordinates(:,1),coordinates(:,2),coordinates(:,3),xyz(:,1),xyz(:,2),xyz(:,3),color)
        end
    end
end
    