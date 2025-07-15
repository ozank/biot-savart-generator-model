%---------------------------------------------------
%  NAME:      Plot Axial Race Track Winding.m
%  WHAT:      Plots the conventional race track winding and calculates the magnetic field using the Biot-Savart Calculation Method
%  AUTHOR:    Ozan Keysan (06/2025)
%  REQUIRED:  BSmag Toolbox 20150407 (for magnetic field calculations)
%----------------------------------------------------


% coils are placed symmetrically around Z-axis
Coil_z_min = -0.5 *coil_to_coil_gap * machine.Nstacks; %minimum coil z(-) for the outer stack


for s = 1:machine.Nstacks+1 %Create machine.Nstacks+1 number of coils in the axial direction, midplane Z=0

winding_coordinates_stack = winding_coordinates;

winding_coordinates_stack(:,3) = Coil_z_min + (s-1)*coil_to_coil_gap; %z coordinate of -z stack

% 2D Rotation about a point
% https://academo.org/demos/rotation-about-point/

for n = 1:Npoles_radial %rotate the modules in radial direction by rotating the original coil

winding_coordinates_rotated = winding_coordinates_stack;
%rotation angle
rotation_angle = 1 * (n-1) * machine.pole_angle;
winding_coordinates_rotated(:,1) = winding_coordinates(:,1)*cosd(rotation_angle) - winding_coordinates(:,2)*sind(rotation_angle);
winding_coordinates_rotated(:,2) = winding_coordinates(:,2)*cosd(rotation_angle) + winding_coordinates(:,1)*sind(rotation_angle);

%add rotated coils ( on the + z axis)

% Change direction of current for alternating poles!!!
if mod(n, 2) == 0
  % pole is even -> Reverse current polarity
[BSmag] = BSmag_add_filament(BSmag,winding_coordinates_rotated,-I,dGamma2); %add reversed coil pole

else
  % pole is odd -> Keep current polarity positive

[BSmag] = BSmag_add_filament(BSmag,winding_coordinates_rotated,I,dGamma2); % add normal polarity pole

end

end

end

%End Winding Coil Loops

% for n = 1:( Npoles_radial/2)  % Add end windings on one side with half of the pole number per module * module number
% 
% %rotation angle per end winding over Z-axis
% rotation_angle = (2*(n-1) + 0.5) * machine.pole_angle;
% 
% %Rotate initial end winding coordinates
% end_winding_coordinates_rotated = end_winding_coordinates;
% 
% %rotation angle
% end_winding_coordinates_rotated(:,1) = end_winding_coordinates(:,1)*cosd(rotation_angle) - end_winding_coordinates(:,2)*sind(rotation_angle);
% end_winding_coordinates_rotated(:,2) = end_winding_coordinates(:,2)*cosd(rotation_angle) + end_winding_coordinates(:,1)*sind(rotation_angle);
% 
% %Add End Winding Coil
% [BSmag] = BSmag_add_filament(BSmag,end_winding_coordinates_rotated,I,dGamma2);
% 
% %add the other end coil (mirror over the Z-axis)
% end_winding_coordinates_rotated(:,3) = -1 * end_winding_coordinates_rotated(:,3);     % Take symmetry over X-Y plane, by multlipying Z coordinate by -1
% [BSmag] = BSmag_add_filament(BSmag,end_winding_coordinates_rotated,I,dGamma2);
% 
% end

