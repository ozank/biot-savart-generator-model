%---------------------------------------------------
%  NAME:      Calculate Electrical Parameters .m
%  WHAT:     Calculates the electrical parameters starting from a single turn coil
%  then populates is to induced voltage and phase current per phase using the machine parameters
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

% Prerequsite: Call the main function and flux calculation before this function


%% Induced Voltage Calculations
stator.induced_voltage_per_coil = stator.N_turns * (2 * pi()/sqrt(2)) * machine.f_electrical * flux;  %[Vrms] d Flux_peak/dt, converted to RMS value 

stator.induced_voltage_per_phase = stator.induced_voltage_per_coil * stator.Nseries;  %[Vrms], Induced voltage per stator phase (per stack)



%% Phase Current Calculations

stator.winding_area = stator.coil_thickness * stator.coil_width * 10^6; %[mm^2], stator winding area, cross section
stator.conductor_area = stator.winding_area * stator.fill_factor / stator.N_turns; %[mm^2], stator conductor (per turn area)
stator.current_per_coil = stator.current_density * stator.conductor_area; %[Arms], current in the conductor of single coil

stator.phase_current = stator.current_per_coil * stator.Nparalel; %[A,rms], phase current per phase (per axial stack) of machine