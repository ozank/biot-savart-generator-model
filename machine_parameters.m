%---------------------------------------------------
%  NAME:      Machine Parameters.m
%  WHAT:      Main Machine Parameters for the Axial Flux Machine
%  AUTHOR:    Ozan Keysan (06/2025)
%----------------------------------------------------
% 
%   Inputs:
%       No inputs
% 
%   Outputs:
%      Independent and Dependent Parameters for the machine (please see the
%      HTSHAAM generator documentation for parameter details)
%
%   Notes:
%       - Dependent and Independent Parameters to be added for radius and
%       other parameters


%% Main Machine Parameters

%machine.Npole = 120;                        % Number of Poles
machine.pole_angle = 360 / machine.Npole;   % One pole pitch angle in degrees

machine.Nrpm = 10;                          %Rated rotational speed of the generator (RPM)
machine.f_electrical = (machine.Nrpm /60)* (machine.Npole/2);   %Machine induced voltage frequency (Hz)

%Number of axial stacks
%machine.Nstacks = 3; % Number of generator stator stacks in the axial direction

stator.Ncoil = machine.Npole * 3/4;        %Number of stator coils per stage, assumes 3/4 relation
machine.Nphase = 3;                         %Number of phases, default to 3
stator.Ncoil_per_phase = stator.Ncoil / machine.Nphase;  %Number of stator coils per phase (per axial stack)

%% Mean (Radius) Parameters

% Rotor
% Mean radius (independent)
%HTS.R_mean = 3000/1000;     %meter, Mean radius of the HTS rotor
%HTS.coil_length = 300/1000;  %meters, Straight section length of rotor HTS coil
%Outer and other parameters are calculated below using number of turns,
%wire dimensions etc.

%Stator
stator.R_mean = HTS.R_mean;  %can be modified but by default, equal to rotor mean radius
%stator.coil_length = 300/1000;  %meters, Straight section length of the stator coils, independent but correlated with HTS.coil_length, stator.R_inner, Stator.R_outer
stator.coil_length = HTS.coil_length; %make stator length equal to HTS length (temporary)

%% HTS Winding Type
%Only use one type by commenting out the other one

HTS.winding_type = 'race_track'; % HTS rotor winding is conventional trapezoidal race track winding
%HTS.winding_type = 'wave';   %HTS rotor winding is in the shape of wave winding (i.e. potato masher)


%% Rotor HTS Parameters and Excitation Current
HTS.current = 550; %A   HTS Coil Current (DC), per tape
%HTS.N_turns = 100; %    Number of turns of HTS coil (per layer)
HTS.N_layers = 2;  %    Number of layers, default is 2 for double pancake configurations

%Define the total current for biot savart conductor
I = HTS.current * HTS.N_turns * HTS.N_layers; %[A] filament current in the biot savart model

%% Stator Coil Parameters

%stator.N_turns = 50;             %Number of turns in a single stator coil
stator.Nparalel = 1;            %Number of paralel connected stator winding coils per phase (default to 1)
stator.Nseries = stator.Ncoil_per_phase/stator.Nparalel;  %Number of series connected stator winding per phase 

% Current Density & Fill Factor
%stator.current_density = 7;    %[A/mm^2, RMS], Current density in stator windings (calculated in the conductor area)
stator.fill_factor = 0.8;       %Fill factor (independent or calculated according to stator coil cross section area)

%Litz Wire Parameters
stator.Litz_N_strands = 50;   %Number of strands of the Litz wire, used for eddy/proximity loss calculations, a value of 50-200 seems enough for DD wind turbines
stator.Litz_fill_factor = 1; % Ignored for now, the overall fill factor is included in the stator.fill_factor, can be separated later on


%% Stator Coil Width and Height Calculations
%stator.coil_width_to_coil_pitch_ratio = 0.4;  %Coil width to pole pitch ratio (at mean radius), <0.5 but considering inner radius 0.45 is a more realistic limit

stator.coil_angle = 360 / stator.Ncoil;   % One stator coil pole pitch angle in degrees
stator.coil_pitch = stator.R_mean * (stator.coil_angle *pi() /180); %[m], Coil pitch at mean radius

%Stator Coil Width & Thickness
stator.coil_width = stator.coil_pitch * stator.coil_width_to_coil_pitch_ratio ;    %meters, Width of the stator coil on one side
%stator.coil_thickness = 25/1000;    %meters, Thickness of single stator coil

stator.R_bending = 30/1000;     %meters, Bending radius of the stator coils at the four outer corners

%% Stator Outer Inner Dimensions

stator.R_outer = stator.R_mean + 0.5*stator.coil_length + stator.coil_width;  %Mean diameter plus half of straight section plus coil width
stator.R_inner = stator.R_mean - 0.5*stator.coil_length - stator.coil_width;  %Mean diameter minus half of straight section minus coil width

%% HTS Tape Dimensions

HTS.tape_width = 12; %[mm], single HTS tape width
HTS.tape_thickness = 0.2; %[mm] single HTS tape thickness

%% HTS Coil Parameters, Rotor Outer Inner Dimensions
%Secondary parameters
HTS.pole_pitch = 2 * pi * HTS.R_mean / machine.Npole;  %[m] Magnetic pole pitch of the rotor HTS coils
HTS.distance_to_next = 10/1000; %[m], distance between two adjacent HTS coils at mean radius
HTS.gap_pancake = 1/1000; %[m], Gap between two pancake coils of HTS, default to 1 mm

HTS.coil_pitch = HTS.pole_pitch - HTS.distance_to_next; %[m] Coil pitch at mean radius of HTS coils
HTS.coil_width = HTS.tape_thickness * HTS.N_turns / 1000; %[m] Coil Width of the HTS coil 
%Geometry check
% HTS.coil_width < 0.5 *HTS.coil_pitch !!

HTS.coil_thickness = (HTS.tape_width/1000 * HTS.N_layers) + HTS.gap_pancake; %[m] coil thickness (axial direction of HTS coil (double pancake))

%will need to be defined with respect to other parameters (to be changed)
HTS.R_outer = HTS.R_mean + 0.5 * HTS.coil_length + HTS.coil_width; %[m], HTS Outer radius, similar calculation to stator
HTS.R_inner = HTS.R_mean - 0.5 * HTS.coil_length - HTS.coil_width; %[m], HTS Inner radius, similar calculation to stator


%% Air Gap, Axial Gap, Cryostat Thickness
machine.airgap_mechanical = 4.5/1000; %[m], mechanical airgap between stator and rotor
machine.airgap_cryostat = 10/1000; %[m], Gap required for cryostat, includes vacuum thickness in axial length, and the thicknes of the cryostat in axial direction (if exists): HTS surface to mechanical gap starting point
machine.airgap_magnetic = machine.airgap_cryostat + machine.airgap_mechanical; %[m] Magnetic airgap (i.e. HTS surface to stator coil surface, depends on cryostat dimensions and mechanical airgap)

machine.airgap_HTS_to_HTS = 2 * machine.airgap_magnetic + stator.coil_thickness; %[m], Axial distance between to HTS coil surfaces between stacks (axial direction)

coil_to_coil_gap = machine.airgap_HTS_to_HTS + HTS.coil_thickness; %[m], HTS filament to filament distance in the axial direction (between stacks)
%coil_to_coil_gap = 0.1; %manual Z-direction ignoring coil thicknes(ie filament), mid point of the coils should be considered


%% Race Track Coil Parameters
%Coil Angle (ratio of coil width to pole pitch at the mean diameters
% can be be connected to other machine parameters and name can be modified
%TO BE REMOVED!!, moved to coordinate calculation function
HTS.coil_angle = (HTS.coil_pitch / HTS.pole_pitch) * machine.pole_angle; %[degrees], coil span angle considering the gaps between adjacent coils

%% Wave Winding Shape Adjustments
HTS.pole_outer_ratio = 1; %pole ratio of the outer pole edge, default to 1  0 < value < 2
HTS.pole_inner_ratio = 1; %pole ratio of the outer pole edge, default to 1  0 < value < 2

machine.Npole_per_module = 8; %Number of poles per cryostat module, default to 8, should be an even number

%Only for the wave winding type machine, radius of the wave winding bottom
%return path 
%TO BE ADJUSTED!
HTS.R_bottom = HTS.R_inner - 0.3;   %meter, Bottom radius of the HTS rotor for wave winding (to be adjusted)


%% End Winding Coordinates
% To be modified Refer to Figure
% Rotation angle of the end winding coils
end_winding_type = 3; %End winding type: 0= no end winding, 1 = only vertical, 3 = vertical with diagonal ones (three in total)
end_winding_rotation_angle= 45;

end_winding_gap = 0.1;   %Z distance from the last stage of coils (the name can be modified)
end_winding_Z = 0.5 *coil_to_coil_gap * machine.Nstacks + end_winding_gap ;  %Z-coordinate of the end winding
