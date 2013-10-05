classdef Test_Plotting < matlab.unittest.TestCase
    methods (Test)
        function test_H8(testCase)
            dim = 2;
            a = 2; m = dim;
            b = 1; n = dim;
            c = 3; p = dim;
            E = 2;                                % Youngs modulus
            nu = 0.3;                               % Poisson's ratio
            rho = 7840;                           % Density
            sides = [a,b,c];
            mesh = Mesh.meshgen('Mech_H8',sides,[m n p],E,nu,rho);
            node_function = (1:mesh.nnodes)';
            figure
            mesh.plot_node_function(node_function)
        end
        function test_Q4(testCase)
            dim = 2;
            a = 2; m = dim;
            b = 1; n = dim;
            E = 2;                                % Youngs modulus
            nu = 0;                               % Poisson's ratio
            rho = 7840;                           % Density
            sides = [a,b];
            mesh = Mesh.meshgen2D('Mech_Q4',sides,[m n],E,nu,rho);
            node_function = (1:mesh.nnodes)';
            figure
            mesh.plot_node_function(node_function)
        end
    end
end