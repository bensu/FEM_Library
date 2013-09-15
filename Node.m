classdef Node < hgsetget
    properties 
        id
        coordinates
        parent
    end
    methods
        function obj = Node(id,coordinates,parent)
            set(obj,'id',id);
            set(obj,'coordinates',coordinates);
            set(obj,'parent', parent);
        end
        function dofsid = ldofsid(obj)
            %returns the list with it's dof 
            mesh = obj.get('parent');
            ndofpernode = mesh.ndofspernode();
            nodeids = obj.get('id');
            dofsid = (1+ndofpernode*(nodeids-1):ndofpernode*nodeids);
            if length(dofsid)~=ndofpernode
                error('DEBUG NODE.GDOFSID')
            end
        end
            
    end
end