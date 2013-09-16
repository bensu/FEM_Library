classdef Test_Mech_Q4 < matlab.unittest.TestCase
  
    methods (Test)
        function patch_test(test_case)
            for dim = 1:2
                for coordnum = 1:2
                    plot_on = false;
                    %% Mesh & Material
                    % Cubes dimension a,b,c
                    % Number of elements in each cube edge.
                    a = 2; m = dim;
                    b = 1; n = dim;
                    
                    E = 2;                                % Youngs modulus
                    nu = 0;                               % Poisson's ratio
                    rho = 7840;                           % Density
                    sides = [a,b];
                    mesh = Mesh.meshgen2D('Mech_Q4',sides,[m n],E,nu,rho);
                    
                    % Break Mesh Symmetry
                    delta = 0.1;
                    mesh.randominner(delta);    % Moves an inner node by delta
                    
                    %% BC
                    valuebc = 0;
                    face = Face.new_face(mesh,coordnum,valuebc);
                    bc = BC(mesh.sdof());
                    bc.simplysupportedface(face,coordnum,mesh);
                    
                    %% LOADS
                    q = 2;
                    coordload = coordnum;
                    valueload = sides(coordload);
                    vector = zeros(2,1);
                    vector(coordload) = q;
                    face = Face.new_face(mesh,coordload,valueload);
                    loads = Loads(mesh.sdof());
                    loads.qinface(mesh,face,vector);
                    
                    %% CASE
                    patch = FemCase(mesh,mesh.sdof(),mesh.sdof());
                    patch.set('loads',loads);
                    patch.set('bc',bc);
                    patch.solve();
                    patch.getMaxStressEle()
                    SA = patch.get('StressArrayEle');
                    %subplot(3,1,coordnum)
                    if plot_on
                        patch.plot()
                    end
                    
                    %% Check
                    
                    dis = patch.get('displacements');
                    maxval = dis.maxvalue();
                    tol = 1e-12;
                    sigma_theorical = q*sides(coordnum)/E;
                    test_case.verifyEqual(abs(maxval-sigma_theorical)<tol,true)
                    
                end
            end 
            patch.plot([])
        end
    end
end
