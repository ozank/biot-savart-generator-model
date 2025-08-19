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

%% Initialize Biot Savart Model
% Clear workspace
clear all, close all, clc
BSmag = BSmag_init(); % Initialize BSmag analysis

%BIOT SAVART MODEL SETTINGS
%Discrete steps along coil during, biot savart calculations, has a direct
%effect on calculation time
dGamma2 = 1e-2; % filament max discretization step [m], default to 1 cm  

%% Get Material Properties
material_constants; % Load material constants

%% Optimization Parameters

% Number of Inputs
number_of_inputs = 2; %Number of optimization inputs

%Inputs
%inputs = [3 1];
%Number of poles
%Current Density

%Lower and Upper Bounds for the optimization Inputs
bounds = [90 150;     %Number of Poles
          2  10]      %J (current density) A/mm^2

% Lower Bounds for optimization inputs
lower_bounds = bounds(:,1)   % First column assigned to lower bounds

%Upper Bounds for optimization inputs
upper_bounds = bounds(:,2)   % Second column assigned to upper bounds


%% Flux Per Pole (with Biot Savart Model)
%Determine number of data points for airgap flux density calculations
data_point_angle= 20;  % number of data points in the tangential directions (through angle)
data_point_radius = 50; %number of data points in the radial (radius) direction


rng default
nvars = 2;
opts = optimoptions(@gamultiobj,'PlotFcn','gaplotpareto');
[xga,fvalga,~,gaoutput] = gamultiobj(@(x)optimization_cost(x),nvars,[],[],[],[],lower_bounds,upper_bounds,[],opts);