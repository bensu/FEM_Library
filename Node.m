classdef Node < hgsetget
    properties 
        id
        coordinates
        v3
    end
    methods
        function obj = Node(id,coordinates)
            set(obj,'id',id);
            set(obj,'coordinates',coordinates);
        end
    end
end