classdef Q4 < Element
    methods 
        function obj = Q4(nodes,material,element_dofs_in, ...
                                element_displacements_in)
            require(length(nodes)==4,'Needs 8 nodes')
            require(length(nodes(1).get('coordinates'))==3, ...
                'Mesh should be 3D');
            obj = obj@Element(nodes,material,element_dofs_in, ...
                                element_displacements_in);
        end  
        %% Q4 functions
        % should be inherited from Q4 element
        function [shapeQ4,dNdxiQ4_out,dNdetaQ4_out] = shapefunctions(element,xi,eta)

            %------------------------------------------------------------------------
            %  Purpose:
            %     compute isoparametric four-node Quadilateral shape functions
            %     and their derivatves at the selected (integration) point
            %     in terms of the natural coordinate
            %
            %  Snetanopsis:
            %     [shapeQ4,dhdrQ4,dhdsQ4]=shapefunctions(etai,eta)
            %
            %  Variable Description:
            %     shapeQ4 - shape functions for four-node element
            %     dhdrQ4 - derivatives of the shape functions w.r.t. r
            %     dhdsQ4 - derivatives of the shape functions w.r.t. s
            %     etai - r coordinate value of the selected point
            %     eta - s coordinate value of the selected point
            %
            %  Notes:
            %     1st node at (-1,-1), 2nd node at (1,-1)
            %     3rd node at (1,1), 4th node at (-1,1)
            %------------------------------------------------------------------------
            
            shapeQ4 = element.shapeQ4(xi,neta);
            dNdxiQ4_out = element.dNdetaQ4(xi,eta);
            dNdetaQ4_out = element.dNdnetaQ4(xi,eta); 
        end 
        function shapeQ4_out = N(element,xi,eta)
            %  Notes:
            %     1st node at (-1,-1), 3rd node at (-1,1) 
            %     4th node at (1,1), 2nd node at (1,-1)
            require(isscalar(xi)&isscalar(eta),'Both xi and eta should be scalar')
            % shape functions
            shapeQ4_out = zeros(1,4);
            shapeQ4_out(1) = 0.25*(1-xi)*(1-eta);
            shapeQ4_out(3) = 0.25*(1+xi)*(1-eta);
            shapeQ4_out(2) = 0.25*(1-xi)*(1+eta);
            shapeQ4_out(4) = 0.25*(1+xi)*(1+eta);
        end        
        function dNdxiQ4_out = dN_dxi(element,xi,eta)
            % derivatives
            require(isscalar(xi)&isscalar(eta),'Both xi and eta should be scalar')
            dNdxiQ4_out = zeros(1,4);
            dNdxiQ4_out(1) = -0.25*(1-eta);
            dNdxiQ4_out(2) = 0.25*(1-eta);
            dNdxiQ4_out(3) = -0.25*(1+eta);
            dNdxiQ4_out(4) = 0.25*(1+eta);
        end
        function dNdetaQ4_out = dN_deta(element,xi,eta)
            % derivatives
            require(isscalar(xi)&isscalar(eta),'Both xi and eta should be scalar')
            dNdetaQ4_out = zeros(1,4);
            dNdetaQ4_out(1) = -0.25*(1-xi);
            dNdetaQ4_out(2) = -0.25*(1+xi);
            dNdetaQ4_out(3) = 0.25*(1-xi);
            dNdetaQ4_out(4) = 0.25*(1+xi);
        end
        function dN_out = DN(element,xi,eta)
            dN_out = 
    end
end
    