classdef Test_Mech_ShellQ4 < matlab.unittest.TestCase
    
    methods (Test)
        function patch_test(test_case)
            plot_on = false;
            a = 1; b = 1;
            thickness = 0.1;
            mesh_patch_test = Mesh.ShellQ4Mesh(a,b,thickness,2,2,1,1,0.3,1);
            coords = mesh_patch_test.get('coordinates');
            coords(5,1:2) = coords(5,1:2) + [0.1, -0.1];        %Symmetry
            mesh_patch_test.set('coordinates',coords);
            patch = FemCase(mesh_patch_test);
            bc = patch.get('bc');
            bc.fix_node_by_id(1); %Punto empotrado
            bc.support_node_by_id(2,2);
            bc.support_node_by_id(3,2);
            loads = patch.get('loads');
            sigma = 2;
            F = sigma*a*thickness/4;
            loads.get('node_component').edit_component_by_id(7,[0 F 0 0 0]);
            loads.get('node_component').edit_component_by_id(8,[0 2*F 0 0 0]);
            loads.get('node_component').edit_component_by_id(9,[0 F 0 0 0]);
            
            %% CASE
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