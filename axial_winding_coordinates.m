%---------------------------------------------------
%  NAME:      Axıal Winding Coordinates.m
%  WHAT:      Defines the initial coordinates and variables for the
%  conventional axial winding (trapezoidal) HTSAAM machine
%  AUTHOR:    Ozan Keysan (06/2025)
%----------------------------------------------------
% 
%   Inputs:
%       No inputs
% 
%   Outputs:
%       Race track coil parameters and coordinate references required (please see the generator documentation for parameter details
%

% Axial Adjustment of Rotor Coils
coil_Z_offset = 0; % coil is placed on X-Y plane, coil Z offset defines Z coordinate


%CREATE TRAPEZOIDAL WINDING DIMENSIONS

%refer to the figure
%create a zero matrix with the correct order
winding_coordinates = zeros (5,3);  %Initialize default matrix for coordinates of the wave winding  

   
%Point A
%Inner radius of the first pole, first coil segment
angle =  (machine.pole_angle - coil_angle)/2;
R = HTS.R_inner;

winding_coordinates(1,:) =  [   R * cosd(angle),
                                R * sind(angle),
                                coil_Z_offset];

%Point B
%Outer radius of the first pole, first coil segment
angle =  (machine.pole_angle - coil_angle)/2;
R = HTS.R_outer;

winding_coordinates(2,:) =  [   R * cosd(angle),
                                R * sind(angle),
                                coil_Z_offset];

%Point C
%outer radius, first pole, second coil segment
angle =  (machine.pole_angle - coil_angle)/2 + coil_angle;
R = HTS.R_outer;

winding_coordinates(3,:) =  [   R * cosd(angle),
                                R * sind(angle),
                                coil_Z_offset];

%Point D
%inner radius, first pole, second coil segment
angle =  (machine.pole_angle - coil_angle)/2 + coil_angle;
R = HTS.R_inner;

winding_coordinates(4,:) =  [   R * cosd(angle),
                                R * sind(angle),
                                coil_Z_offset];

%Point A
% Finish at the starting point
% Inner radius of the first pole, first coil segment
angle =  (machine.pole_angle - coil_angle)/2;
R = HTS.R_inner;

winding_coordinates(5,:) =  [   R * cosd(angle),
                                R * sind(angle),
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
angle = coil_angle;
R = HTS.R_outer;

end_winding_coordinates(3,:) =  [   R * cosd(angle),
                                    0,
                                    end_winding_Z + R * sind(angle)];


%pointD
%Inner radius, outer edge of the end winding
angle = coil_angle;
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