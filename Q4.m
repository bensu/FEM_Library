classdef Q4 < Element
    properties (Constant)
        nnodes = 4;
        node_connectivity = [1 2;1 3;2 4;3 4];
        face_connectivity = [1 2 4 3];
        node_local_coords = [-1 -1;-1 1;1 -1;1 1];
    end
    methods 
        function ele_out = Q4(id_in,n_node_dofs,n_element_dofs,nodes_in,material_in)
            require(length(nodes_in(1).get('coordinates'))==2, ...
                'Mesh should be 2D');
            ele_out = ele_out@Element(id_in,n_node_dofs,n_element_dofs, ...
                                nodes_in,material_in);
        end  
        %% Shape Functions
        % should be inherited from Q4 element
        function [shapeQ4,dNdxiQ4_out,dNdetaQ4_out] = shapefunctions(element,local_coords)

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
            [xi, eta] = element.components(local_coords);
            shapeQ4 = element.shapeQ4(xi,eta);
            dNdxiQ4_out = element.dNdetaQ4(xi,eta);
            dNdetaQ4_out = element.dNdnetaQ4(xi,eta); 
        end 
        function shapeQ4_out = N(element,local_coords)
            [xi, eta] = element.components(local_coords);
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
        function dNdxiQ4_out = dN_dxi(element,local_coords)
            [xi, eta] = element.components(local_coords);
            % derivatives
            require(isscalar(xi)&isscalar(eta),'Both xi and eta should be scalar')
            dNdxiQ4_out = zeros(1,4);
            dNdxiQ4_out(1) = -0.25*(1-eta);
            dNdxiQ4_out(2) = 0.25*(1-eta);
            dNdxiQ4_out(3) = -0.25*(1+eta);
            dNdxiQ4_out(4) = 0.25*(1+eta);
        end
        function dNdetaQ4_out = dN_deta(element,local_coords)
            [xi, eta] = element.components(local_coords);
            % derivatives
            require(isscalar(xi)&isscalar(eta),'Both xi and eta should be scalar')
            dNdetaQ4_out = zeros(1,4);
            dNdetaQ4_out(1) = -0.25*(1-xi);
            dNdetaQ4_out(2) = -0.25*(1+xi);
            dNdetaQ4_out(3) = 0.25*(1-xi);
            dNdetaQ4_out(4) = 0.25*(1+xi);
        end
        function dN_dmu(element,local_coords)
            error('dN_dmu should not be caled on a 2D element')
        end
        function dN_out = DN(element,local_coords)
            dN_out = [element.dN_dxi(local_coords);element.dN_deta(local_coords)];
        end
        
    end
end
    