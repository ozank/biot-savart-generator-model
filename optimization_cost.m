%---------------------------------------------------
%  NAME:      Optimization Cost.m
%  WHAT:      Cost funtion for the genetic algoritm:
%  - Magnetostatic analysis (With Biot Savart Model)
%  - Electrical Equivalent Circuit Parameters
%  - Loss Calculations
%  - Mass/Cost Calculation
%
%  Inputs: Defined by the Pre Optimization Function
%  Outputs: Cost function (suitable for multi objective optimization)
%
%  AUTHOR:    Ozan Keysan (08/2025)
%  REQUIRED:  BSmag Toolbox 20150407 (for magnetic field calculations)
%----------------------------------------------------

function cost = optimization_cost(inputs)

BSmag = BSmag_init(); % Initialize BSmag analysis


%BIOT SAVART MODEL SETTINGS
%Discrete steps along coil during, biot savart calculations, has a direct
%effect on calculation time
dGamma2 = 1e-2; % filament max discretization step [m], default to 1 cm  

%% Get Material Properties
material_constants; % Load material constants


%Input arrangements

machine.Npole = 4*floor(inputs(1)/4);   %Number of Poles (make sure it is a multiple of 4
stator.current_density = inputs(2);     % Current Density (A/mm^2)


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

% Default Optimization Function Tries to Minimize
cost = zeros(1,2); % allocate output depending on number of parameters

cost(1) =  - machine.efficiency;  % [0 1] Efficiency of the machine
cost(2) =  - machine.P_output;    % [W] Electrical Output Power

%Close all figures
close all

