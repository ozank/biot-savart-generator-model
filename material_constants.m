%---------------------------------------------------
%  NAME:      Material Constants.m
%  WHAT:      Material constants and properties to be used for other functions
%  AUTHOR:    Ozan Keysan (07/2025)
%----------------------------------------------------
% 
%   Inputs:
%       No inputs
% 
%   Outputs:
%      Material Properties

%% Physical Contstants
material.mu0 = 4*pi()*10^-7;  % [A/m], mu_0, permeability of free space

material.copper_resistivity_at_room_temperature =  1.68 * 10^-8; %[Ohm/m] copper resistivity at room temperature (20 C)
material.aluminium_resistivity_at_room_temperature =  2.65 * 10^-8; %[Ohm/m] aluminium resistivity at room temperature (20 C)

%https://cirris.com/temperature-coefficient-of-copper/
material.copper_temperature_constant = 0.00393;  %Temperature coefficient for copper
material.aluminium_temperature_constant = 0.004308;  %Temperature coefficient for aluminium


%% Operating Temperature
stator.operating_temperature = 80;      %[C], Assumed operating temperature of the stator windings, can be linked to a thermal model later on


%% Effective Resistivity at Operating Temperature
%Stator Winding Material : Copper
%material.winding_resistivity = (1 + material.copper_temperature_constant * (stator.operating_temperature - 20)) ...
%                                * material.copper_resistivity_at_room_temperature; %[Ohm/m], effective resistivity at operating temperature 

%Hypercooled aluminium resistivity
material.winding_resistivity =  2.12 * 10^-10;  %hypercooled aluminium at 35 K


%Stator Winding Material: Aluminium (uncomment the below core if the material is aluminium
%material.winding_resistivity = (1 + material.aluminium_temperature_constant * (stator.operating_temperature - 20))...
%                                * material.aluminium_resistivity_at_room_temperature; %[Ohm/m], effective resistivity at operating temperature 


%% Mass & Cost Properties

%% Density of materials
material.copper_density = 8960; %[kg/m^3] Density of copper
material.aluminium_density = 2700; %[kg/m^3] Density of aluminium
material.HTS_density = 8000; %[kg/m^3] Density of HTS tape

material.winding_density = material.copper_density; %[kg/m^3], chose between aluminium and copper windings depending on the design

%% Cost of Materials (to be adjusted)
material.copper_cost = 35; % [$/kg], unit cost of copper wire
material.aluminium_cost = 8; % [$/kg], unit cost of aluminium wire

material.winding_cost = material.copper_cost; %[GBP/kg],  chose between aluminium and copper windings depending on the design

% HTS cost per meter (from Shuangrong 08/2025)
%  4mm  40-50 $/kg
% 12 mm 120-150 $/kg from Super power
% a reduction is assumed for bulk orders
material.HTS_cost = 50; %[$/m], or [$/kg], to be adjusted, cost of HTS YBCO tape

