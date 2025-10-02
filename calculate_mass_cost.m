%---------------------------------------------------
%  NAME:      Calculate Mass & Cost.m
%  WHAT:     Calculates the active material mass and cost of the generator 
%  AUTHOR:    Ozan Keysan (08/2025)
%----------------------------------------------------
% 
%   Inputs:
%       Machine Dimensions Flux per pole and machine dimensions
%       Current Density
% 
%   Outputs:
%       Mass and Cost of the following active materials
% - Copper Wire
% - HTS Tape

% Prerequsite: Call the main function and electrrical parameter calculations before this function

% TO DO: Structural Mass Components to be added

%% Calculate Winding Mass
% length of the copper, multiplied by area by copper density
% Overall mass (including axial stacks)
stator.mass = material.winding_density *...
              (stator.mean_turn_length * stator.conductor_area * 10^-6 * stator.N_turns) * ...
              stator.Ncoil * machine.Nstacks; %[kg], overall mass of the winding (copper/aluminium) for all coils


%% HTS Length and Mass Calculations

%Different Calculations for the Race Track and Wave Winding Types

%Filament length for mean HTS (per coil)(for mass/cost calculations)
% For race track coils or the end winding type coils
HTS.mean_turn_length = 2*(filament_HTS.R_outer - filament_HTS.R_inner) + 2 * (HTS.coil_pitch - HTS.coil_width); %[m] mean turn length of single HTS coil


% For race track coils
if strcmp(HTS.winding_type, 'race_track')  %Draw the race track winding

    % Total number of HTS Coils
    %End windings type can be (0, 1 or 3 as stated in the machine parameters)
    
    %configuration below uses double pancake (HTS.N_layers)
    % A double pancake coil is equals to single HTS coil count
    HTS.Ncoils_total =  machine.Npole * (machine.Nstacks + 1) + 2 * (machine.Npole/2) * end_winding_type;
    
    HTS.length_total = HTS.N_turns * HTS.N_layers * HTS.mean_turn_length * HTS.Ncoils_total; %[m] Total length of HTS cable that is required for the design


else
%% HTS Usage Calculations for the Wave Winding

%Length of the loop (1 turn) in one wave winding module

HTS.module_turn_length = (filament_HTS.R_inner - filament_HTS.R_bottom) * 2 ... % From points A to B (including on each side)
                          + sqrt((winding_coordinates(3,1) - winding_coordinates(4,1))^2 + (winding_coordinates(3,2) - winding_coordinates(4,2))^2 ) * machine.Npole_per_module ... %Points B and C of the wave winding (vertical connections), calculate from the Hypotenuse length
                          + HTS.pole_outer_ratio * machine.pole_angle * (pi/180) * filament_HTS.R_outer * (machine.Npole_per_module/2) ...  %outer edge on the Router, between points C and D
                          + HTS.pole_inner_ratio * machine.pole_angle * (pi/180) * filament_HTS.R_inner * (machine.Npole_per_module/2 - 1) ...  %innerr edge on the Rinner, between points E and F
                          + (machine.pole_angle * (machine.Npole_per_module - 1)) * (pi/180) * filament_HTS.R_bottom; % Bottom connection from module end to point A 

%% HTS usage for wave winding
HTS.length_total =  HTS.N_turns * HTS.N_layers * HTS.module_turn_length * (machine.Npole / machine.Npole_per_module) * (machine.Nstacks + 1)  ... %HTS used in wave winding modules
                  + HTS.N_turns * HTS.N_layers * HTS.mean_turn_length * 2 * (machine.Npole/2) * end_winding_type; % HTS used in end windings

end


% Calculate HTS Mass
HTS.mass = material.HTS_density * HTS.length_total * (HTS.tape_width /1000) * (HTS.tape_thickness/1000); %[kg], total mass of the HTS material


%% Material Cost Calculations

% Stator Copper/Aluminium Material Cost Calculation
stator.cost = stator.mass *  material.winding_cost; %[$] material cost of stator windings

%HTS Material Cost
% Can be calculated either per kg or per meter
HTS.cost = HTS.length_total * material.HTS_cost; %[$] material cost of the HTS winding

