classdef Mat < hgsetget
    properties 
        id
        type
        C
        density
    end
    methods
        function obj = Mat(id,type,E,nu,rho)
            set(obj,'type',type);
            if strcmp(type,'Iso')
                C = E/((1+nu)*(1-2*nu))*[1-nu nu   nu   0          0          0          
                       nu   1-nu nu   0          0          0  
                       nu   nu   1-nu 0          0          0
                       0    0    0    (1-2*nu)/2 0          0
                       0    0    0    0          (1-2*nu)/2 0
                       0    0    0    0          0          (1-2*nu)/2];
            end
            set(obj,'C',C);
            set(obj,'id',id);
            set(obj,'density',rho);
        end

    end
end

    