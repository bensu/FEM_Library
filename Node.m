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
            ndofs = mesh.get('dofs_per_node');
            nodeids = obj.get('id');
            dofsid = (1+ndofs*(nodeids-1):ndofs*nodeids);
            if length(dofsid)~=ndofs
                error('DEBUG NODE.GDOFSID')
            end
        end
            
    end
end