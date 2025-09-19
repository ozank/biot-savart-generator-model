%---------------------------------------------------
%  NAME:      Optimization Cost With Penalty.m
%  WHAT:      Cost funtion for the genetic algoritm:
%  - Magnetostatic analysis (With Biot Savart Model)
%  - Electrical Equivalent Circuit Parameters
%  - Loss Calculations
%  - Mass/Cost Calculation
%
%  Inputs: Defined by the Pre Optimization Function
%  Penalty: Penalty functions added for efficiency etc.
%  Outputs: Single Cost function (suitable for singleobjective optimization)
%
%  AUTHOR:    Ozan Keysan (08/2025)
%  REQUIRED:  BSmag Toolbox 20150407 (for magnetic field calculations)
%----------------------------------------------------

function cost = optimization_cost_with_penalty(inputs)

BSmag = BSmag_init(); % Initialize BSmag analysis


%BIOT SAVART MODEL SETTINGS
%Discrete steps along coil during, biot savart calculations, has a direct
%effect on calculation time
dGamma2 = 1e-2; % filament max discretization step [m], default to 1 cm  

%% Get Material Properties
material_constants; % Load material constants


%Input arrangements
%Make sure they are aligned with the input upper and lower limits in GA
%optimization

machine.Npole = 4*inputs(1);   %Number of Poles divided by 4  (make sure it is a multiple of 4
%machine.Npole = 4*floor(inputs(1)/4);   %Number of Poles (make sure it is a multiple of 4
stator.current_density = inputs(2);     % Current Density (A/mm^2)
HTS.coil_length = inputs(3);
HTS.N_turns = inputs(4);
stator.N_turns = inputs(5);
stator.coil_width_to_coil_pitch_ratio = inputs(6);
stator.coil_thickness = inputs(7);
machine.Nstacks = inputs(8);


%% Get Machine Parameters
% Don't forget to comment out any input parameters.
machine_parameters;      %load machine parameters

%% RACE TRACK COIL WINDING
%Get Race track axial machine Winding Coordinates 

axial_winding_coordinates;
Npoles_radial = 8; % Number of modules to be simulated in the radial direction, default 8
 
% Add windings for the axial race track winding
plot_axial_race_track_winding;

%% Electrical Machine Parameter Estimations

%% Flux Per Pole (with Biot Savart Model)
%Determine number of data points for airgap flux density calculations
data_point_angle= 20;  % number of data points in the tangential directions (through angle)
data_point_radius = 50; %number of data points in the radial (radius) direction

calculate_flux_per_pole;  % Outputs flux per pole and maximum B values

%% Calculate Electrical Parameters
% Get induced voltage, current, resistance etc
calculate_electrical_parameters; 

%% Calculate Power Losses and Efficiency
% Get power output, losses(conduction, eddy), efficiency
calculate_power_efficiency;

%% Calculate mass cost function
% active material cost, material prices
calculate_mass_cost;


%% COST FUNCTION ADJUSTMENTS

%% Penalty Functions

%Limits for Penalty
limit.efficiency = 0.95;   % Minimum required efficiency
limit.P_output = 20e6;      % [W], Required power output

%reset penalties
penalty.efficiency = 0; 
penalty.P_output = 0; 
penalty.total = 0; 

%Efficiency Limit Penalty
if (machine.efficiency < limit.efficiency)
    penalty.efficiency = 1e10 * (penalty.efficiency - machine.efficiency)^2;
else
    %penalty.efficiency = 0; 
end    

%Power Output Limit Penalty
if (machine.P_output < limit.P_output)
    penalty.P_output = 1 * (penalty.P_output - machine.P_output)^2;
else
    %penalty.P_output = 0; 
end    


%% Overall Penalty Values

penalty.total = penalty.efficiency + penalty.P_output;

%Objective Function

cost = HTS.cost + stator.cost + penalty.total;

%Close all figures
close all

