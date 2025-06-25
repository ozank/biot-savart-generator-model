%---------------------------------------------------
%  NAME:      Machine Parameters.m
%  WHAT:      Main Machine Parameters for the axial flux machine
%  AUTHOR:    Ozan Keysan (06/2025)
%----------------------------------------------------
% 
%   Inputs:
%       No inputs
% 
%   Outputs:
%      Independent and Dependent Parameters for the machine (please see the
%      generator documentation for parameter details)
%
%   Notes:
%       - Dependent and Independent Parameters to be added for radius and
%       other parameters

%Clear previous parameters
%clear all;


%Independent Machine Parameters

machine.Npole = 120;                        % Number of Poles
machine.pole_angle = 360 / machine.Npole;   % One pole pitch angle in degrees

% Main Radius Definitions
HTS.R_mean = 2525/1000;     %meter, Mean radius of the HTS rotor

HTS.R_outer = 2690/1000;    %meter, Outer radius of the HTS rotor (dependent, to be adjusted)
HTS.R_inner = 2360/1000;    %meter, Inner radius of the HTS rotor (dependent, to be adjusted)
HTS.R_bottom = 2175/1000;   %meter, Bottom radius of the HTS rotor (dependent, to be adjusted)

%HTS Parameters and Excitation Current
HTS.current = 225; %A   HTS Coil Current (DC), per tape
HTS.N_turns = 970; %    Number of turns of HTS coil (per layer)
HTS.N_layers = 2;  %    Number of layers, default is 2







%Wave Winding Shape Adjustments
HTS.pole_outer_ratio = 1; %pole ratio of the outer pole edge, default to 1  0 < value < 2
HTS.pole_inner_ratio = 1; %pole ratio of the outer pole edge, default to 1  0 < value < 2


machine.Npole_per_module = 8; %Number of poles per cryostat module, default to 8, should be an even number


