classdef Test_Mech_Q4 < matlab.unittest.TestCase
  
    methods (Test)
        function patch_test(test_case)
            for dim = 2
                for coordnum = 2
                    plot_on = true;
                    %% Mesh & Material

                    % Plate dimension a,b,thickness
                    % Number of elements in each edge.
                    a = 2; m = dim;
                    b = 1; n = dim;
                    thickness = 1;
                    
                    E = 2;                                % Youngs modulus
                    nu = 0;                               % Poisson's ratio
                    rho = 7840;                           % Density
                    sides = [a,b];
                    mesh = Mesh.meshgen2D('Mech_Q4',sides,[m n],E,nu,rho);
                    % Break Mesh Symmetry
                    delta = 0.1;
                    mesh.random_inner(delta);    % Moves an inner node by delta
                    
                    %% BC
                    valuebc = 0;
                    face = Face.new_face(mesh,coordnum,valuebc);
                    bc = BC(mesh.nnodes,mesh.get('dofs_per_node'), ...
                            mesh.nnel, mesh.get('dofs_per_ele'));
                    bc.simply_supported_face(face,coordnum,mesh);
                    
                    
                    %% LOADS
                    
                    
                    
                    loads = Loads(mesh.nnodes,mesh.get('dofs_per_node'), ...
                            mesh.nnel, mesh.get('dofs_per_ele'));
%                     loads.qinface(mesh,face,vector);
%                     loads.node_function

%                   %% NOT GENERALIZED FOR D
                    sigma = 2;
                    F = sigma*thickness*sides(mod(coordnum,2)+1)/4;
                    q1 = zeros(2,1); q1(coordnum) = F;                    
                    loads.get('node_component').edit_component_by_id(7,q1);
                    loads.get('node_component').edit_component_by_id(8,2*q1);
                    loads.get('node_component').edit_component_by_id(9,q1)
                                    
                    %% CASE
                    patch = FemCase(mesh,bc,loads);
                    patch.solve();
                    stress_array = patch.node_stress_array;
                    %subplot(3,1,coordnum)
                    if plot_on
                        figure
                        patch.plot_displacement(coordnum)
                    end
                    
                    %% Check
                    
                    max_stress = max(max(stress_array));
                    tol = 1e-12;
                    test_case.verifyEqual(abs(max_stress-sigma)<tol,true)
                end
            end 
            
        end
    end
end
