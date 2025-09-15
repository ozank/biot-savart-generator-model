%---------------------------------------------------
%  NAME:      Post Process Optimization.m
%  WHAT:      Read from the optimization solutions to run the main function
% to calculate main machine parameters
%
%  AUTHOR:    Ozan Keysan (09/2025)
%----------------------------------------------------

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


%% MAIN Function Call

main;