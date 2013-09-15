classdef BC < VectorField
    properties

    end

    methods
        function obj = BC(sdof)
            nodelist = true(sdof,1);
            obj = obj@VectorField(nodelist);
        end
        function obj = addnodeconstraint(obj,nodeid,coordnum)
            previous = obj.xyzout();
            previous(nodeid,coordnum) = false;
            obj.set('nodelist',previous);
        end
        function obj = addfixedpoint(obj,nodeid)
            obj.addnodeconstraint(nodeid,1);
            obj.addnodeconstraint(nodeid,2);
            obj.addnodeconstraint(nodeid,3);
        end
        function obj = addfixedface(obj,new_face)
            previous = obj.xyzout();
            nodelist = new_face.get('nodelist');
            for i = 1:length(nodelist)
                previous(nodelist(i),:) = false(1,3);
            end
            obj.set('nodelist',VectorField.xyzin(previous));
        end
        function addfaceconstraint(obj,face,coordnum,mesh)
            nodelistids = mesh.facenodelist(face);
            for i = 1:length(nodelistids)
                obj.addnodeconstraint(nodelistids(i),coordnum);
            end
        end
        function simplysupportedface(obj,face,coordnum,mesh) %missing a point
            coordinates = mesh.get('coordinates');
            facelist = mesh.facenodelist(face);
            %first fixed point
            nodeid = facelist(1);
            obj.addfixedpoint(nodeid);
            %face
            obj.addfaceconstraint(face,coordnum,mesh);
            u = zeros(3,1);
            u(coordnum) = 1;
            iteration = 1;
            log = 0;
            iterationlimit = 1;
            while and(log == 0,iterationlimit<100)
                iteration = iteration + 1;
                coordnum2 = find(cross(u,coordinates(nodeid,:) - coordinates(nodeid+1,:)));
                if ~isempty(coordnum2)
                    obj.addnodeconstraint(nodeid+1,coordnum2);
                    log = 1;
                else nodeid = nodeid + 1;
                end
            end
          
        end  
        function plot(obj,coordinates,color)
            hold on
            xyz = ~obj.xyzout();
            quiver3(coordinates(:,1),coordinates(:,2),coordinates(:,3),xyz(:,1),xyz(:,2),xyz(:,3),color)
        end
    end
end
    