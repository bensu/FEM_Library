classdef Mech_Shell < Mechanics
    methods
        function C_out = C(element)
            C_out = element.C@Mechanics;
            C_out(3,:) = 0;
            C_out(:,3) = 0;
        end
        function Ndevs_out = DNsparse(element,local_coords)
            [~, ~, mu] = element.components(local_coords);
            dim = element.dim;
            ndofs = element.get('dofs_per_node');
            t_at_node = element.t_at_node;
            dN = [element.dN_dxi(local_coords); element.dN_deta(local_coords)];
            N = element.N(local_coords);
            u_matrix = element.u_matrix;
            total_nodes = element.nnodes;
            Ndevs_out = zeros(dim^2,total_nodes*ndofs);  % MAGIC 5
            for node = 1:total_nodes
                aux = zeros(dim^2,ndofs);
                for j = 1:dim
                    aux2 = zeros(dim,ndofs);
                    aux2(1:2,j) = dN(:,node);
                    aux2(1:2,4:5) = 0.5*t_at_node(node)*mu*dN(:,node)*u_matrix(j,:,node);
                    aux2(3,4:5) = 0.5*t_at_node(node)*N(node)*u_matrix(j,:,node);
                    aux(index_range(dim,j),:) = aux2;
                end
                Ndevs_out(:,index_range(ndofs,node)) = aux;
            end
        end
    end
end