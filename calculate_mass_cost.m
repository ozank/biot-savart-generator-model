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
% For race track coils

% Total number of HTS Coils
%End windings are calculated using 3 end coils (2 angled + 1  vertical),
%needs to be adjusted for different configuration
%configuration below uses double pancake (HTS.N_layers)
HTS.Ncoils_total =  machine.Npole * HTS.N_layers * machine.Nstacks + 2 * (machine.Npole/2) * 3;

HTS.length_total = HTS.N_turns *HTS.mean_turn_length * HTS.Ncoils_total; %[m] Total length of HTS cable that is required for the design

HTS.mass = material.HTS_density * HTS.length_total * (HTS.tape_width /1000) * (HTS.tape_thickness/1000); %[kg], total mass of the HTS material

%Todo
%Wave winding to be added

%% Material Cost Calculations

% Stator Copper/Aluminium Material Cost Calculation
stator.cost = stator.mass *  material.winding_cost; %[GBP] material cost of stator windings

%HTS Material Cost
% Can be calculated either per kg or per meter

HTS.cost = HTS.mass * material.HTS_cost; %[GBP] material cost of the HTS winding

