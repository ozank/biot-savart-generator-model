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
%  AUTHOR:    Ozan Keysan (08/2025), modified (09/2025)
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
HTS.R_mean = inputs(9);
%Extra parameters for the wave winding
HTS.distance_inner_to_bottom = inputs(10);
HTS.pole_outer_ratio = inputs(11);
HTS.pole_inner_ratio = inputs(12);

%% Get Machine Parameters
% Don't forget to comment out any input parameters.
machine_parameters;      %load machine parameters

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
    Npoles_radial = 8; % Number of modules to be simulated in the radial direction, default 8
 
    % Add windings for the axial race track winding
    plot_axial_race_track_winding;

    %% Electrical Machine Parameter Estimations
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


%% COST FUNCTION ADJUSTMENTS

%% Penalty Functions

%Limits for Penalty
limit.efficiency = 0.99;   % Minimum required efficiency
limit.P_output = 2e6;      % [W], Required power output
%Extra Limits
limit.Vphase_max = 1200;  %[Vrms] per phase voltage

%reset penalties
penalty.efficiency = 0; 
penalty.P_output = 0; 
penalty.V_phase = 0;
penalty.total = 0; 

%Efficiency Limit Penalty
if (machine.efficiency < limit.efficiency)
    penalty.efficiency = 1e10 * (limit.efficiency - machine.efficiency)^2;
else
    %penalty.efficiency = 0; 
end    

%Power Output Limit Penalty
if (machine.P_output < limit.P_output)
    penalty.P_output = 1 * (limit.P_output - machine.P_output)^2;
else
    %penalty.P_output = 0; 
end    

%Phase Voltage Limit
if (stator.induced_voltage_per_phase > limit.Vphase_max)
    penalty.V_phase = 1e2 * (limit.Vphase_max - stator.induced_voltage_per_phase)^2;
else
    %penalty.V_phase = 0; 
end    


%% Overall Penalty Values

penalty.total = penalty.efficiency + penalty.P_output + penalty.V_phase;

%Objective Function
%Minimum Material Cost
cost = HTS.cost + stator.cost + penalty.total;

%Minimum Active Material Mass
%cost = HTS.mass + stator.mass + penalty.total;


%Close all figures
close all

