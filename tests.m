clc
clear all
clear clases

%% Runs the test suit for the whole Library
%  Every Test should be added sequentially.

addpath('./tests')
% 
% patch_mech_H8 = Test_Mech_H8;
% run(patch_mech_H8)
% 
% patch_mech_Q4 = Test_Mech_Q4;
% run(patch_mech_Q4)
% 
% patch_mech_ShellQ4 = Test_Mech_ShellQ4 ;
% run(patch_mech_ShellQ4)

patch_ME_ShellQ4 = Test_ME_ShellQ4 ;
run(patch_ME_ShellQ4)

% plot_tests = Test_Plotting;
% run(plot_tests)