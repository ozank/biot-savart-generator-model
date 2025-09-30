%---------------------------------------------------
%  NAME:      Main.m
%  WHAT:      Main function for Generator Design that completes the following:
%  - Magnetostatic analysis (With Biot Savart Model)
%  - Electrical Equivalent Circuit Parameters
%  - Loss Calculations
%  - Mass/Cost Calculation
%
%  AUTHOR:    Ozan Keysan (08/2025)
%  REQUIRED:  BSmag Toolbox 20150407 (for magnetic field calculations)
%----------------------------------------------------

%% Initialize Biot Savart Model
%clear all, close all, clc
BSmag = BSmag_init(); % Initialize BSmag analysis

%BIOT SAVART MODEL SETTINGS
%Discrete steps along coil during, biot savart calculations, has a direct
%effect on calculation time
dGamma2 = 1e-2; % filament max discretization step [m], default to 1 cm  

% %Optimization Outputs
% machine.Npole = 240 % 4*solution(1);   %Number of Poles divided by 4  (make sure it is a multiple of 4
% %machine.Npole = 4*floor(inputs(1)/4);   %Number of Poles (make sure it is a multiple of 4
% stator.current_density = 6 %solution(2);     % Current Density (A/mm^2)


%% Get Machine Parameters
machine_parameters;      %load machine parameters
%small_machine_parameters;      %load small machine parameters

%Input arrangements


%% Get Material Properties
material_constants; % Load material constants

%Determine number of data points for airgap flux density calculations
%Has a direct effect on computation time and accuracy
data_point_angle= 15;  % number of data points in the tangential directions (through angle)
data_point_radius = 40; %number of data points in the radial (radius) direction

%% Get Winding coordinates

% It is possible to draw two types of windings: Wave winding and
% conventional race track coils, Please comment out the unwanted type, and
% use ONLY one of the winding types

if strcmp(HTS.winding_type, 'race_track')  %Draw the race track winding
    %% RACE TRACK COIL WINDING
    %Get Race track axial machine Winding Coordinates 

     axial_winding_coordinates;
    Npoles_radial = 6; % Number of modules to be simulated in the radial direction, default 8
 
    % Add windings for the axial race track winding
    plot_axial_race_track_winding;

    %% Flux Per Pole 
    %Calculate Flux Per Pole
    calculate_flux_per_pole;  % Outputs flux per pole and maximum B values

else         %Draw wave winding
    %% WAVE WINDING
    %Get Wave Winding Coordinates
    
    wave_winding_coordinates;
    Nmodules_radial = 2; % Number of modules to be simulated in the radial direction, default 3

    % Add windings for the wave winding
    plot_wave_winding;

    %% Flux Per Pole 
    %Calculate Flux Per Pole
    calculate_flux_per_pole_wave_winding; %Flux per pole calculations for the wave winding

end

%% Calculate Electrical Parameters
% Get induced voltage, current, resistance etc
calculate_electrical_parameters; 

%% Calculate Power Losses and Efficiency
% Get power output, losses(conduction, eddy), efficiency
calculate_power_efficiency;

%% Calculate mass cost function
% active material cost, material prices
calculate_mass_cost;