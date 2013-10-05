classdef (Abstract) Physics < hgsetget
    
%     properties (Abstract, Constant)
%         node_dofs_id    % Dofs Ids this Physical Property takes per node
%         ele_dofs_id     % Dofs Ids this Physical Property takes per element
%     end
%     properties (Dependent)
%         n_node_dofs     % Total number of Dofs this Physical Property takes per node
%         n_ele_dofs      % Total number of Dofs this Physical Property takes per ele
%     end


    methods (Abstract)
        B_out = B(element,xi,eta,mu);
        C_out = C(element,xi,eta,mu);
        K_out = K(element,gauss_order);
    end
    
end