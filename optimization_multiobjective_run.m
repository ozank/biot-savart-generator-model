%---------------------------------------------------
%  NAME:      Optimization Run.m
%  WHAT:      Optimization code and the codes that prepares the variables for optimization.

% Single objective Genetic Algorithm 
%  Cleans workspace, defines variables required by the optimization such as
% - upper and lower bounds
% - optimization options
% - number of variables etc. 

%  AUTHOR:    Ozan Keysan (10/2025)
%----------------------------------------------------

%% Clear Workspace
clear all, close all, clc

%% Plot Options
% Disabling plots can make the optimization a lot faster
% It is recommended to make this variable zero for the optimization loop

global plot_figures
plot_figures = 0; % If 0 do not plot figures, if = 1 then plot figures 

%% Initialize Biot Savart Model
BSmag = BSmag_init(); % Initialize BSmag analysis

%BIOT SAVART MODEL SETTINGS
%Discrete steps along coil during, biot savart calculations, has a direct
%effect on calculation time
dGamma2 = 1e-2; % filament max discretization step [m], default to 1 cm  

%% Flux Per Pole (for the Biot Savart Model)
%Determine number of data points for airgap flux density calculations
data_point_angle= 15;  % number of data points in the tangential directions (through angle)
data_point_radius = 40; %number of data points in the radial (radius) direction

%% Get Material Properties
material_constants; % Load material constants

%% Optimization Parameters

% Number of Inputs
number_of_inputs = 12; %Number of optimization inputs

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
% 9 - HTS.R_mean

% Extra Parameters for the wave winding type
%10 - HTS.distance_inner_to_bottom
%11 - HTS.pole_outer_ratio;
%12 - HTS.pole_inner_ratio;

%Lower and Upper Bounds for the optimization Inputs
bounds = [4 4;      % Number of Poles/4 ( the value divided by 4) due to simulation constraints
          10 40;       % J (current density) A/mm^2
          0.075 0.2;    % 3- HTS.coil_length
          10 120;     % 4- HTS.N_turns
          10 200;     % 5- stator.N_turns
          0.2 0.5;   % 6- stator.coil_width_to_coil_pitch_ratio
          0.015 0.03;      % 7- stator.coil_thickness    
          1  1;          %8- machine.Nstacks
          0.25  0.25; %];          % 9- HTS.R_mean   
%Extra parameters for the wave winding type
          0.2 0.5;   % 10- HTS.distance_inner_to_bottom
          0.8 1.2;   % 11- HTS.pole_outer_ratio;
          0.8 1.2];  %12 - HTS.pole_inner_ratio;


% Lower Bounds for optimization inputs
lower_bounds = bounds(:,1);   % First column assigned to lower bounds

%Upper Bounds for optimization inputs
upper_bounds = bounds(:,2);   % Second column assigned to upper bounds

%IntegerConditions
% Define the parameters that take integer value ( for shorter converge
% time)
IntCon = [1, 4, 5, 8]; % The first variable is integer or [1,3] first and third inputs are integers


%% Optimization Settings
% Set nondefault solver options
options = optimoptions("gamultiobj","PopulationSize",250, ...
                            "MaxGenerations",100,...
                            "MaxStallGenerations",20,... 
                            "UseParallel",true, ...
                            "ConstraintTolerance",0.1, ...
                            "FunctionTolerance",0.1,...
                            "Display","iter", ...
                            "PlotFcn",["gaplotpareto","gaplotscores","gaplotstopping","gaplotspread","gaplotrankhist","gaplotrankhist"]);


%% RUN Genetic Algorithm
% Solve
% Main optimization function: Optimization_multi_cost
[solution,objectiveValue] = gamultiobj(@optimization_multi_cost,number_of_inputs,...
    [],[],[],[],lower_bounds,upper_bounds,[],IntCon,options);

% Clear variables
clearvars options

%Plot Pareto Front
figure()
plot(objectiveValue(:,2),-objectiveValue(:,1),...
     'o','MarkerSize',8','LineWidth',2)
ylabel('Efficiency (%)')
xlabel('Active Material Mass (kg)')
xlabel('Required HTS Tape Length (m)')
fontsize(12,"points")
grid on