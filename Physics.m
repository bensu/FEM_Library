classdef (Abstract) Physics < hgsetget

    methods (Abstract)
        B_out = B(element,xi,eta,mu);
        C_out = C(element,xi,eta,mu);
        K_out = K(element,gauss_order);
    end
end