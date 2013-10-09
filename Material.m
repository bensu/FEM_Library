classdef Material < hgsetget
    properties 
        id
        type
        Young_Modulus           % [Pa]
        Poisson_Coefficient     % []
        Density                 % [kg/m3]
        Permitivity3
        Piezo13
        Piezo23        
    end
    methods
        function mat = Material(id,type,Young_Modulus_in, ...
                Poisson_Coefficient_in, Density_in,varargin)
%           mat = Material(id,type,E,nu,rho)
%           mat = Material(id,type,E,nu,rho,permitivity3,piezo13,piezo23)
            set(mat,'id',id);
            set(mat,'type',type);
            set(mat,'Young_Modulus',Young_Modulus_in);
            set(mat,'Poisson_Coefficient',Poisson_Coefficient_in);
            set(mat,'Density',Density_in);
            if length(varargin) >= 3 
                set(mat,'Permitivity3',varargin{1})
                set(mat,'Piezo13',varargin{2})
                set(mat,'Piezo23',varargin{3})
            end
        end

    end
end

    