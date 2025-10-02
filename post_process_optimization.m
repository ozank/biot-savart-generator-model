%---------------------------------------------------
%  NAME:      Post Process Optimization.m
%  WHAT:      Read from the optimization solutions to run the main function
% to calculate main machine parameters
%
%  AUTHOR:    Ozan Keysan (09/2025)
%----------------------------------------------------

%% Plot Options
% Disabling plots can make the optimization a lot faster
% It is recommended to make this variable zero for optimization loop
% but it can be set to 1 to generate post processing figures

global plot_figures
plot_figures = 1; % If 0 do not plot figures, if = 1 then plot figures 

%% Variables
%Read solution variable from optimization outputs
%Input arrangements
%Make sure they are aligned with the actual optimization cost function

machine.Npole = 4*solution(1);   %Number of Poles divided by 4  (make sure it is a multiple of 4
%machine.Npole = 4*floor(inputs(1)/4);   %Number of Poles (make sure it is a multiple of 4
stator.current_density = solution(2);     % Current Density (A/mm^2)
HTS.coil_length = solution(3);
HTS.N_turns = solution(4);
stator.N_turns = solution(5);
stator.coil_width_to_coil_pitch_ratio = solution(6);
stator.coil_thickness = solution(7);
machine.Nstacks = solution(8);
HTS.R_mean = solution(9);
%extra parameters for the wave winding
%can be commented out for the race track winding
HTS.distance_inner_to_bottom = solution(10);
HTS.pole_outer_ratio = solution(11);
HTS.pole_inner_ratio = solution(12);

%% MAIN Function Call

main;