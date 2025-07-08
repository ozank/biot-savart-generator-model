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

machine.Nrpm = 10;                          %Rated rotational speed of the generator (RPM)
machine.f_electrical = (machine.Nrpm /60)* (machine.Npole/2);   %Machine induced voltage frequency (Hz)

%Number of axial stacks
machine.Nstacks = 3; % Number of generator stator stacks in the axial direction


machine.Ncoil = machine.Npole * 3/4;        %Number of stator coils per stage, assumes 3/4 relation
machine.Nphase = 3;                         %Number of phases, default to 3


coil_to_coil_gap = 0.1; %ignoring coil thicknes, mid plane of the coils should be considered

% Main Radius Definitions
HTS.R_mean = 2525/1000;     %meter, Mean radius of the HTS rotor

%will be defined with respect to other parameters
HTS.R_outer = 2690/1000;    %meter, Outer radius of the HTS rotor (dependent, to be adjusted)
HTS.R_inner = 2360/1000;    %meter, Inner radius of the HTS rotor (dependent, to be adjusted)
HTS.R_bottom = 2000/1000;   %meter, Bottom radius of the HTS rotor (dependent, to be adjusted)

%HTS Parameters and Excitation Current
HTS.current = 225; %A   HTS Coil Current (DC), per tape
HTS.N_turns = 970; %    Number of turns of HTS coil (per layer)
HTS.N_layers = 2;  %    Number of layers, default is 2


%Wave Winding Shape Adjustments
HTS.pole_outer_ratio = 1; %pole ratio of the outer pole edge, default to 1  0 < value < 2
HTS.pole_inner_ratio = 1; %pole ratio of the outer pole edge, default to 1  0 < value < 2


machine.Npole_per_module = 8; %Number of poles per cryostat module, default to 8, should be an even number

%Stator Coil Parameters

stator.N_turns = 1;             %Number of turns in a single stator coil
stator.coil_length = 300/1000;  %meters, Straight section length of the stator coils, independent but correlated with HTS.coil_length, stator.R_inner, Stator.R_outer
stator.Nparalel = 1;            %Number of paralel connected stator winding coils per phase
stator.coil_thickness = 30/1000;%meters, Thickness of single stator coil
stator.fill_factor = 0.6;       %Fill factor (independent or calculated according to stator coil cross section area)

stator.coil_width = 30/1000;    %meters, Width of the stator coil on one side
stator.R_bending = 30/1000;     %meters, Bending radius of the stator coils at the four corners

%Should be calculated according to dimensions
stator.R_outer = HTS.R_outer;  %to be changed
stator.R_inner = HTS.R_inner;  %to be changed
stator.R_mean = HTS.R_mean;  %to be changed
