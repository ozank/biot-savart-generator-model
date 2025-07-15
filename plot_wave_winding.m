%---------------------------------------------------
%  NAME:      Plot Wave Winding.m
%  WHAT:      Plots the wave windings and calculates the magnetic field using the Biot-Savart Calculation Method
%  AUTHOR:    Ozan Keysan (06/2025)
%  REQUIRED:  BSmag Toolbox 20150407 (for magnetic field calculations)
%----------------------------------------------------


% coils are placed symmetrically around Z-axis
Coil_z_min = -0.5 *coil_to_coil_gap * machine.Nstacks; %minimum coil z(-) for the outer stack

%Create stack coils 
% coils are placed symmetrically around Z-axis

for s = 1:machine.Nstacks+1 %Create machine.Nstacks+1 number of coils in the axial direction, midplane Z=0

winding_coordinates_stack = winding_coordinates;

winding_coordinates_stack(:,3) = Coil_z_min + (s-1)*coil_to_coil_gap; %z coordinate of -z stack

% 2D Rotation about a point
% https://academo.org/demos/rotation-about-point/

for n = 1:Nmodules_radial %rotate the modules in radial direction by rotating the original coil

winding_coordinates_rotated = winding_coordinates_stack;
%rotation angle
rotation_angle = 1 * (n-1) * machine.pole_angle * machine.Npole_per_module;
winding_coordinates_rotated(:,1) = winding_coordinates(:,1)*cosd(rotation_angle) - winding_coordinates(:,2)*sind(rotation_angle);
winding_coordinates_rotated(:,2) = winding_coordinates(:,2)*cosd(rotation_angle) + winding_coordinates(:,1)*sind(rotation_angle);

%add rotated coils ( on the + z axis)
[BSmag] = BSmag_add_filament(BSmag,winding_coordinates_rotated,I,dGamma2);

end

end

%End Winding Coil Loops

for n = 1:(machine.Npole_per_module/2 * Nmodules_radial)  % Add end windings on one side with half of the pole number per module * module number

%rotation angle per end winding over Z-axis
rotation_angle = (2*(n-1) + 0.5) * machine.pole_angle;

%Rotate initial end winding coordinates
end_winding_coordinates_rotated = end_winding_coordinates;

%rotation angle
end_winding_coordinates_rotated(:,1) = end_winding_coordinates(:,1)*cosd(rotation_angle) - end_winding_coordinates(:,2)*sind(rotation_angle);
end_winding_coordinates_rotated(:,2) = end_winding_coordinates(:,2)*cosd(rotation_angle) + end_winding_coordinates(:,1)*sind(rotation_angle);

%Add End Winding Coil
[BSmag] = BSmag_add_filament(BSmag,end_winding_coordinates_rotated,I,dGamma2);

%add the other end coil (mirror over the Z-axis)
end_winding_coordinates_rotated(:,3) = -1 * end_winding_coordinates_rotated(:,3);     % Take symmetry over X-Y plane, by multlipying Z coordinate by -1
[BSmag] = BSmag_add_filament(BSmag,end_winding_coordinates_rotated,I,dGamma2);

end
