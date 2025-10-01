%---------------------------------------------------
%  NAME:      Pre Optimization.m
%  WHAT:      Codes that prepares the variables for optimization toolbox.
%  Cleans workspace, defines variables required by the optimization such as
% - upper and lower bounds
% - optimization options
% - number of variables etc. 

%  AUTHOR:    Ozan Keysan (08/2025)
%  REQUIRED:  BSmag Toolbox 20150407 (for magnetic field calculations)
%----------------------------------------------------

%% Clear Workspace
clear all, close all, clc

%% Plot Options
% Disabling plots can make the optimization a lot faster
% It is recommended to make this variable zero for optimization loop

%global plot_figures

%plot_figures = 0; % If 0 do not plot figures, if = 1 then plot figures 

%% Initialize Biot Savart Model

BSmag = BSmag_init(); % Initialize BSmag analysis


%BIOT SAVART MODEL SETTINGS
%Discrete steps along coil during, biot savart calculations, has a direct
%effect on calculation time
dGamma2 = 1e-2; % filament max discretization step [m], default to 1 cm  

%% Get Material Properties
material_constants; % Load material constants

%% Optimization Parameters

% Number of Inputs
number_of_inputs = 9; %Number of optimization inputs

%Inputs
%inputs = [3 1];
% 1- machine.Npole  %Number of Poles divided by 4  (make sure it is a multiple of 4
% 2- stator.current_density
% 3- HTS.coil_length
% 4- HTS.N_turns
% 5- stator.N_turns
% 6- stator.coil_width_to_coil_pitch_ratio
% 7- stator.coil_thickness
% 8- machine.Nstacks

%Lower and Upper Bounds for the optimization Inputs
bounds = [20 20;      % Number of Poles ( the value divided by 4) due to simulation constraints
          4 8;       % J (current density) A/mm^2
          0.4 0.8;    % 3- HTS.coil_length
          80 160;     % 4- HTS.N_turns
          80 150;     % 5- stator.N_turns
          0.3 0.45;   % 6- stator.coil_width_to_coil_pitch_ratio
          0.03 0.07;      % 7- stator.coil_thickness    
          5  5;          %8- machine.Nstacks
          2  2];          % HTS.R_mean   

% Lower Bounds for optimization inputs
lower_bounds = bounds(:,1);   % First column assigned to lower bounds

%Upper Bounds for optimization inputs
upper_bounds = bounds(:,2);   % Second column assigned to upper bounds

%IntegerConditions
% Define the parameters that take integer value ( for shorter converge
% time)
IntCon = [1, 4, 5, 8]; % The first variable is integer or [1,3] first and third inputs are integers

%% Flux Per Pole (with Biot Savart Model)
%Determine number of data points for airgap flux density calculations
data_point_angle= 20;  % number of data points in the tangential directions (through angle)
data_point_radius = 50; %number of data points in the radial (radius) direction

