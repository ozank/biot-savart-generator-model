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

stator.coil_resistance = stator.N_turns * stator.mean_turn_length * material.winding_resistivity / (stator.conductor_area *10^-6); %[Ohm], DC resistance of single stator coil

stator.phase_resistance = (stator.coil_resistance * stator.Nseries ) / stator.Nparalel; %[Ohm], DC phase resistance of the stator (per stack)



%% Stator Winding Inductance Calculation

% to be added !!
% either analytical approximation, or biot-savart model can be used for
% inductance calculation.



%% Conduction (DC) Losses
% calculates the power loss (I^2*Rdc) for one phase
% can be moved to external function
stator.P_conduction_loss_per_stage = machine.Nphase * stator.phase_current^2 * stator.phase_resistance; %[W] I^2Rdc per stage of the generator

%% Eddy (AC) Losses

% Calculated according to the equations presented in:
% Manolopoulos, C.D., et al.: Litz wire loss performance and optimization for cryogenic windings. 
% IET Electr. Power Appl. 17(4), 487–498 (2023). https://doi.org/10.1049/elp2.12279 

%Also presented in Hongye's paper
%H. Zhang et al., "Optimizing Loss Performance in Aluminum Litz Wires for Cryogenic Electrical Machines," 
% in IEEE Transactions on Energy Conversion, doi: 10.1109/TEC.2025.3545209.

%Ptot = Prms + Pskin + Pprox
%Equation(2) second term is neglected as it is too small if the litz wire
%diameter is smaller than the skin depth

stator.B_max = 2; %[T], maximum flux density acting on the stator coils for eddy loss estimations

%Calculate Skin Depth
material.skin_depth = sqrt (material.winding_resistivity/(pi * machine.f_electrical * material.mu0))*1000; %[mm], skin depth of the stator conductor (at operating temperature), in mm

%Stator skin loss poer stage (quite small and can be neglected for DD wind
%turbines, but can be significant for high speed applications
stator.P_skin_loss_per_stage = stator.P_conduction_loss_per_stage ...
    * ((stator.Litz_N_strands*(stator.Litz_N_strands-1)*(stator.Litz_strand_diameter/stator.Litz_wire_diameter)^2*(stator.Litz_strand_diameter/material.skin_depth)^4)/128); %[W], calculated according to DC losses

stator.P_proximity_loss_per_stage = stator.Ncoil * stator.N_turns * stator.mean_turn_length *...
(stator.Litz_N_strands/8) * (2 * pi * machine.f_electrical)^2 * stator.B_max^2 * (1/ material.winding_resistivity) * pi * (stator.Litz_strand_diameter/2000)^4  ; %[W], proximity loss including all stator coils (in single stage)


stator.P_eddy_loss_per_stage = stator.P_skin_loss_per_stage + stator.P_proximity_loss_per_stage; %[W] AC losses per stage due to proximity and eddy losses


%% Total Power 
machine.power_factor = 1; %Assume the power factor is unity
machine.P_mechanical = machine.Nphase * machine.Nstacks * stator.induced_voltage_per_phase * stator.phase_current * machine.power_factor; %[W], watts total machine mechanical power

machine.P_loss = machine.Nstacks * (stator.P_conduction_loss_per_stage + stator.P_eddy_loss_per_stage); %[W], total losses in the machine
%other losses to be added, cooling budget

machine.P_output = machine.P_mechanical - machine.P_loss; %[W], Net machine total output power in watts

%% Efficiency Calculations
machine.efficiency = machine.P_output / machine.P_mechanical; % Effficiency (0,1)