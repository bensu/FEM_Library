classdef Material < hgsetget
    properties 
        id
        type
        Young_Modulus           % [Pa]
        Poisson_Coefficient     % []
        Density                 % [kg/m3]
    end
    methods
        function obj = Material(id,type,Young_Modulus_in, ...
                Poisson_Coefficient_in, Density_in)
            set(obj,'id',id);
            set(obj,'type',type);
            set(obj,'Young_Modulus',Young_Modulus_in);
            set(obj,'Poisson_Coefficient',Poisson_Coefficient_in);
            set(obj,'Density',Density_in);
        end

    end
end

    