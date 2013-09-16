classdef VectorField < hgsetget
    properties
        dim         % components per node
        nodelist
    end
    methods (Static)
        function nodelist_out = xyzin(nodelist_in)
            %Takes a Nx3 array and returns a (3*N)x1
           nodelist_out = reshape(nodelist_in',[],1);
        end
    end
    methods
        function obj = VectorField(dim_in,nodelist_in)
            set(obj,'nodelist',nodelist_in);
            set(obj,'dim',dim_in);
        end
        function xyz = xyzout(obj)
            %Takes a (3*N)x1 from nodelist and returns a Nx3
            xyz = reshape(obj.get('nodelist')',obj.get('dim'),[])';
        end
        function set.nodelist(obj,nodelist_in)
            if (size(nodelist_in,1)==1 && size(nodelist_in,2)~=3)
                aux = nodelist_in';
            elseif size(nodelist_in,2)==1
                aux = nodelist_in;
            elseif size(nodelist_in,2)==3
                aux = reshape(nodelist_in',[],1);
            end
            obj.nodelist = aux;
        end
        function M = maxvalue(obj)
            M = max(obj.get('nodelist'));
        end
        function plotVF(obj,coordinates,colour)
            hold on
            Vect3 = obj.xyzout();
            quiver3(coordinates(:,1),coordinates(:,2),coordinates(:,3),Vect3(:,1),Vect3(:,2),Vect3(:,3),colour)
        end
    end
end
    