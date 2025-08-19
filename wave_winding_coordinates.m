%---------------------------------------------------
%  NAME:      Wave Winding Coordinates.m
%  WHAT:      Defines the initial coordinates and variables for the wave winding type axial flux machine
%  AUTHOR:    Ozan Keysan (06/2025)
%----------------------------------------------------
% 
%   Inputs:
%       No inputs
% 
%   Outputs:
%       Parameters and coordinate references required (please see the generator documentation for parameter details
%
%   Notes:
%       - Normal race tracks coil coordinates can be added

% Axial Adjustment of Rotor Coils
coil_Z_offset = 0; % coil is placed on X-Y plane, coil Z offset defines Z coordinates

%% TO DO List

% - Adjust the router and Rinner according to filament dimensions
% - Calculate the mean turn lenght for the single module
% Be careful the code below calculates mean length per coil not length per
% module!
%Filament length for mean HTS (per coil)(for mass/cost calculations)
%HTS.mean_turn_length = 2*(filament_HTS.R_outer - filament_HTS.R_inner) + 2 * (HTS.coil_pitch - HTS.coil_width); %[m] mean turn length of single HTS coil



%CREATE WAVE WINDING DIMENSIONS

%refer to the figure
%create a zero matrix with the correct order
winding_coordinates = zeros (((machine.Npole_per_module/2)*4+3),3);  %Initialize default matrix for coordinates of the wave winding  

%Point A
%Inner starting point (at bottom radius)
winding_coordinates(1,:) = [    HTS.R_bottom*cosd(machine.pole_angle*0.5),
                                HTS.R_bottom*sind(machine.pole_angle*0.5),
                                coil_Z_offset];

%Create poles (depending on the number of poles per cryostat)
for p = 0:(machine.Npole_per_module/2)-1
    
%Point B
%Inner radius of the first pole
angle =  machine.pole_angle*(0.5 + 2*p + 0.5 - HTS.pole_inner_ratio*0.5);
R = HTS.R_inner;

winding_coordinates(2+4*p,:) =  [   R * cosd(angle),
                                    R * sind(angle),
                                    coil_Z_offset];

%Point C
%Outer radius of the first pole
angle =  machine.pole_angle*(0.5 + 2*p + 0.5 - HTS.pole_outer_ratio*0.5);
R = HTS.R_outer;

winding_coordinates(3+4*p,:) =  [   R * cosd(angle),
                                    R * sind(angle),
                                    coil_Z_offset];

%Point D
%first pole end at the outer radius
angle =  machine.pole_angle*(0.5 + 2*p + 0.5 + HTS.pole_outer_ratio*0.5);
R = HTS.R_outer;

winding_coordinates(4+4*p,:) =  [   R * cosd(angle),
                                    R * sind(angle),
                                    coil_Z_offset];


%Point E
%first pole end at the inner radius
angle =  machine.pole_angle*(0.5 + 2*p + 0.5 + HTS.pole_inner_ratio*0.5);
R = HTS.R_inner;

winding_coordinates(5+4*p,:) =  [   R * cosd(angle),
                                    R * sind(angle),
                                    coil_Z_offset];

end

%End of Croyostat point at the other end
winding_coordinates((machine.Npole_per_module/2)*4+2,:) =   [HTS.R_bottom*cosd(machine.pole_angle*(0.5+(machine.Npole_per_module-1))),
                                                            HTS.R_bottom*sind(machine.pole_angle*(0.5+(machine.Npole_per_module-1))),
                                                            coil_Z_offset];
%Point K
%Finish at the starting point = point A
winding_coordinates((machine.Npole_per_module/2)*4+3,:) =    [HTS.R_bottom*cosd(machine.pole_angle*0.5),
                                                            HTS.R_bottom*sind(machine.pole_angle*0.5),
                                                            coil_Z_offset];
%****************************************


%CREATE END WINDING SINGLE COIL COORDINATES

end_winding_coordinates = zeros (5,3);  %Initialize default matrix for coordinates of the end wave winding, 5 points

%pointA
%Inner radius of the end winding
angle = 0;
R = HTS.R_inner;

end_winding_coordinates(1,:) =  [   R,
                                    0,
                                    end_winding_Z];
%pointB
%Outer radius of the end winding
angle = 0;
R = HTS.R_outer;

end_winding_coordinates(2,:) =  [   R,
                                    0,
                                    end_winding_Z];

%pointC
%Outer radius, outer edge of the end winding
angle = machine.pole_angle;
R = HTS.R_outer;

end_winding_coordinates(3,:) =  [   R * cosd(angle),
                                    0,
                                    end_winding_Z + R * sind(angle)];


%pointD
%Inner radius, outer edge of the end winding
angle = machine.pole_angle;
R = HTS.R_inner;

end_winding_coordinates(4,:) =  [   R * cosd(angle),
                                    0,
                                    end_winding_Z + R * sind(angle)];


%pointE = point A
%Inner radius of the end winding
angle = 0;
R = HTS.R_inner;

end_winding_coordinates(5,:) =  [   R,
                                    0,
                                    end_winding_Z];