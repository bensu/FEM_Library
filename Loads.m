classdef Loads < Compound_Function
    properties

    end

    methods
        function loads_out = Loads(total_number_of_nodes_in,dofs_per_node_in, ...
                            total_number_of_elements_in, dofs_per_element_in)
            loads_out = loads_out@Compound_Function(0,total_number_of_nodes_in,dofs_per_node_in, ...
                            total_number_of_elements_in, dofs_per_element_in);    
        end
        function qinface1(obj,mesh,face,qvector)
            qvector= reshape(qvector,[],1);
            [elements, coordlist, values] = face.faceinelement();
            for ele = 1:length(elements)
                coordnum = coordlist(ele);
                value = values(ele);
                localv = zeros(3,1);
                localv(coordnum) = value;
                nn = Element.shfn(localv(1),localv(2),localv(3));
                N = zeros(3,24);
                for j = 1:8
                    N(:,3*(j-1)+1:3*j) = nn(j)*eye(3);
                end
                element = mesh.elementcreate(elements(ele));
                aux = ones(3,1);
                aux(coordnum) = 0;
                logs = find(aux);
                jac = element.jacob(localv(1),localv(2),localv(3));
                detjac = det(jac(logs,logs));
                nodesidlist = element.nodeidlist();
                AUX = reshape(detjac*4*N'*qvector,3,[])';
                for nn = 1:length(nodesidlist);
                    obj.addnode(nodesidlist(nn),AUX(nn,:));
                end
            end
            obj.xyzout();
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
