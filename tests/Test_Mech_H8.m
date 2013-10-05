classdef Test_Mech_H8 < matlab.unittest.TestCase
  
    methods (Test)
        function patch_test(test_case)
            for dim = 2
                for coordnum = 1
                    plot_on = true;
                    %% Mesh & Material
                    n_dim = 3;
                    % Cubes dimension a,b,c
                    % Number of elements in each cube edge.
                    a = 2; m = dim;
                    b = 1; n = dim;
                    c = 3; p = dim;
                    
                    E = 2;                                % Youngs modulus
                    nu = 0.3;                               % Poisson's ratio
                    rho = 7840;                           % Density
                    sides = [a,b,c];
                    mesh = Mesh.meshgen('Mech_H8',sides,[m n p],E,nu,rho);
                    
                    % Break Mesh Symmetry
                    delta = 0.1;
                    mesh.randominner(delta);    % Moves an inner node by delta
                    
                    %% BC
                    valuebc = 0;
                    face = Face.new_face(mesh,coordnum,valuebc);
                    bc = BC(mesh.nnodes,mesh.get('n_node_dofs'), ...
                            mesh.nnel, mesh.get('n_element_dofs'));
                    bc.simply_supported_face(face,coordnum,mesh);
                    
                    %% LOADS
                    sigma = 2;
                    coordload = coordnum;
                    valueload = sides(coordload);
                    vector = zeros(3,1);
                    vector(coordload) = sigma;
                    face = Face.new_face(mesh,coordload,valueload);
                    loads = Loads(mesh.nnodes,mesh.get('n_node_dofs'), ...
                            mesh.nnel, mesh.get('n_element_dofs'));
                    loads.qinface(mesh,face,vector);
                    
                    %% CASE
                    patch = FemCase(mesh,bc,loads);
                    patch.solve();
                    patch.getMaxStressEle();
                    SA = patch.get('StressArrayEle');
                    %subplot(3,1,coordnum)
                    if plot_on
                        patch.plot()
                    end
                    
                    %% Check
                    
                    max_stress = max(max(SA));
                    tol = 1e-12;
                    test_case.verifyEqual(abs(max_stress-sigma)<tol,true)
                    
                end
            end 
        end
    end
end
