%---------------------------------------------------
%  NAME:      Calculate Power & Efficiency .m
%  WHAT:     Calculates the power output, electrical losses and efficiency
%  AUTHOR:    Ozan Keysan (07/2025)
%----------------------------------------------------
% 
%   Inputs:
%       Flux per pole and machine dimensions
%       Current Density
% 
%   Outputs:
%       Induced voltage(RMS) per phase (phase to neutral)
%       Per phase current

% Prerequsite: Call the main function and electrrical parameter calculations before this function


%% Stator Winding Resistance Calculation

%%Mean turn length of coil

%stator.mean_turn_length = 1 ; %%[m], meters calculate from flux per pole
%function, only uncomment for manual calculation

stator.coil_resistance = stator.N_turns * stator.mean_turn_length * material.winding_resistivity / (stator.conductor_area *10^-6); %[Ohm], resistance of single stator coil

stator.phase_resistance = (stator.coil_resistance * stator.Nseries ) / stator.Nparalel; %[Ohm], phase resistance of the stator (per stack)



%% Stator Winding Inductance Calculation

% to be added !!
% either analytical approximation, or biot-savart model can be used for
% inductance calculation.



%% Conduction (DC) Losses
% calculates the power loss (I^2*Rdc) for one phase
% can be moved to external function

stator.P_conduction_loss_per_stage = machine.Nphase * stator.phase_current^2 * stator.phase_resistance; %[W] I^2Rdc per stage of the generator

%% Eddy (AC) Losses
% to be calculated according to B distribution and frequency
%at the moment assumed same as DC conduction losses !!

stator.P_eddy_loss_per_stage = stator.P_conduction_loss_per_stage;   %to be changed

%% Total Power 
machine.power_factor = 1; %Assume the power factor is unity
machine.P_mechanical = machine.Nphase * machine.Nstacks * stator.induced_voltage_per_phase * stator.phase_current * machine.power_factor; %[W], watts total machine mechanical power

machine.P_loss = machine.Nstacks * (stator.P_conduction_loss_per_stage + stator.P_eddy_loss_per_stage); %[W], total losses in the machine
%other losses to be added, cooling budget

machine.P_output = machine.P_mechanical - machine.P_loss; %[W], Net machine total output power in watts

%% Efficiency Calculations
machine.efficiency = machine.P_output / machine.P_mechanical; % Effficiency (0,1)