classdef Loads < Compound_Function
    properties

    end

    methods
        function loads_out = Loads(total_number_of_nodes_in,dofs_per_node_in, ...
                            total_number_of_elements_in, dofs_per_element_in)
            loads_out = loads_out@Compound_Function(0,total_number_of_nodes_in,dofs_per_node_in, ...
                            total_number_of_elements_in, dofs_per_element_in);    
        end
       
        function qinface(loads,mesh,face,qvector)
            qvector = reshape(qvector,1,[]);
            [elelist, facecoord, facevalue] = face.faceinelement();
            NewLoads = zeros(loads.number_of_nodes,loads.dofs_per_node);
            localcoords = zeros(1,3);
            for i = 1:length(elelist)
                new_element = mesh.element_create(elelist(i));
                nodes = new_element.nodeidlist();
                gaussn = 2;
                [gaussp, gaussw] = lgwt(gaussn,-1,1);
                for k1 = 1:gaussn
                    localcoords(1) = gaussp(k1);
                    for k2 = 1:gaussn
                        localcoords(2) = gaussp(k2);
                        for k3 = 1:gaussn
                            localcoords(3) = gaussp(k3);
                            
                            localcoords(facecoord(i)) = facevalue(i);
                            xi = localcoords(1); eta = localcoords(2); mu = localcoords(3);
                            
                            jac =  new_element.jacobian([xi,eta,mu]);
                            index = 1:3;
                            index = index(index~=facecoord(i));
                            jacremainder = jac(index,1:3);
                            V1 = jacremainder(1,:);
                            V2 = jacremainder(2,:);
                            dS = norm(cross(V1,V2));
                            weight = gaussw(k1)*gaussw(k2)*gaussw(k3);
                            NewLoads(nodes',:) = NewLoads(nodes',:) + weight*new_element.N([xi,eta,mu])'*qvector*dS;
                        end
                    end
                end
            end
            loads.add_to_node_function(NewLoads/2); %/2 because there was a double in the loop.
        end
        function plot(obj,coordinates,color)
            hold on
            miniscale = 1e-2;
            xyz = obj.xyzout();
            quiver3(coordinates(:,1)+miniscale,coordinates(:,2),coordinates(:,3),xyz(:,1),xyz(:,2),xyz(:,3),color)
        end
                
            
    end
end
