classdef Test_ME_ShellQ4 < matlab.unittest.TestCase
    
    methods (Test)
        function patch_test(test_case)
            plot_on = false;
            E = 123e9; nu = 0; rho = 1;
            e13 = -5; e23 = 0; e33 = 12.5e-9;
            a = 1; b = 1; thickness = 0.01;
            m = 2; n = 2; p = 1;
            mesh_patch_test = Mesh.ShellQ4Mesh(a,b,thickness,m,n,p,E,nu,rho);
            mesh_patch_test.set('material',Material(1,'Piezo',E,nu,rho,e33,e13,e23));
            mesh_patch_test.set('element_type','ME_ShellQ4');
%             mesh_patch_test = Mesh('ME_ShellQ4',mesh_aux.coords, ...
%                 mesh_aux.connect,Material(1,'Piezo',E,nu,rho,e13,e23,e33));  
            
            coords = mesh_patch_test.get('coordinates');
%             coords(5,1:2) = coords(5,1:2) + [0.1, -0.1];        %Symmetry
            mesh_patch_test.set('coordinates',coords);
            
            patch = FemCase(mesh_patch_test);
            
            %% BC
            bc = patch.get('bc');
            bc.get('element_component');
            bc.fix_node_by_id(1);           %Punto empotrado
            bc.support_node_by_id(2,[2 4 5]);
            bc.support_node_by_id(3,[2 4 5]);
            
            bc.support_element_by_id(1:mesh_patch_test.nnel,2); % Fix all Voltages
            
            %% Loads
            loads = patch.get('loads');
            sigma = 2;
%             F = sigma*a*thickness/2;
%             loads.get('node_component').edit_component_by_id(3,[0 F 0 0 0]);
%             loads.get('node_component').edit_component_by_id(4,[0 F 0 0 0]);
            F = sigma*a*thickness/4;
            loads.get('node_component').edit_component_by_id(7,[0 F 0 0 0]);
            loads.get('node_component').edit_component_by_id(8,[0 2*F 0 0 0]);
            loads.get('node_component').edit_component_by_id(9,[0 F 0 0 0]);
            
            %% CASE
            patch.solve();
            patch.get('displacements').node_function
            patch.get('displacements').element_function
            stress_array = patch.node_stress_array
            stress_array(:,7)
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