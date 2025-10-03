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

%% End Winding Coil Loops

%% Vertical End Winding Coils
% for wnd winding type = 1 or = 3, do not plot if there is no end winding

if end_winding_type > 0 % only plot vertical end winding types for types 1 or 3 


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

else
    %do nothing
end

%% End Winding Coil Loops
% (Diagonal  ones (if used), but angle can be adjusted)

if end_winding_type == 3 % only plot end winding types if the diagonal winding type is selected

for n = 1:(machine.Npole_per_module/2 * Nmodules_radial)  % Add end windings on one side with half of the pole number per module * module number

%% 1st Rotated End coil
%Rotate initial end winding coordinates
end_winding_coordinates_temp = end_winding_coordinates;

%first bring the coil midplane to z=0 plane
z_offset = min(end_winding_coordinates_temp(1:end-1,3)); %Get the average Z-coordinates, except the last one as it is identical with the first one

%move coil to z=0 plane
end_winding_coordinates_temp(:,3) = end_winding_coordinates_temp(:,3) - z_offset;

%rotation angle, rotate around Y axis
end_winding_coordinates_rotated(:,1) = end_winding_coordinates_temp(:,1);
end_winding_coordinates_rotated(:,2) = end_winding_coordinates_temp(:,2)*cosd(end_winding_rotation_angle) + end_winding_coordinates_temp(:,3)*sind(end_winding_rotation_angle);
end_winding_coordinates_rotated(:,3) = end_winding_coordinates_temp(:,3)*cosd(end_winding_rotation_angle) - end_winding_coordinates_temp(:,2)*sind(end_winding_rotation_angle);

%move to z position aligned with the vertical end coils (i.e minimum z of
%the vertical end coil)
end_winding_coordinates_rotated(:,3) = end_winding_coordinates_rotated(:,3) + (min(end_winding_coordinates(:,3))-min(end_winding_coordinates_rotated(:,3)));

%radially rotate so that aligned with the coils in the stacks
rotation_angle2 = (2*n-1) * machine.pole_angle - 0.5*(filament_HTS.coil_angle);

%Rotate initial end winding coordinates
end_winding_coordinates_temp = end_winding_coordinates_rotated;

%rotation angle
end_winding_coordinates_rotated(:,1) = end_winding_coordinates_temp(:,1)*cosd(rotation_angle2) - end_winding_coordinates_temp(:,2)*sind(rotation_angle2);
end_winding_coordinates_rotated(:,2) = end_winding_coordinates_temp(:,2)*cosd(rotation_angle2) + end_winding_coordinates_temp(:,1)*sind(rotation_angle2);
end_winding_coordinates_rotated(:,3) = end_winding_coordinates_temp(:,3);

%Add End Winding Coil
[BSmag] = BSmag_add_filament(BSmag,end_winding_coordinates_rotated,I,dGamma2);

%add the other end coil (mirror over the Z-axis)
end_winding_coordinates_rotated(:,3) = -1 * end_winding_coordinates_rotated(:,3);     % Take symmetry over X-Y plane, by multlipying Z coordinate by -1
[BSmag] = BSmag_add_filament(BSmag,end_winding_coordinates_rotated,I,dGamma2);


%% 2nd Rotated End coil
%Rotate initial end winding coordinates
end_winding_coordinates_temp = end_winding_coordinates;

%first bring the coil midplane to z=0 plane
z_offset = min(end_winding_coordinates_temp(1:end-1,3)); %Get the average Z-coordinates, except the last one as it is identical with the first one

%move coil to z=0 plane
end_winding_coordinates_temp(:,3) = end_winding_coordinates_temp(:,3) - z_offset;

%rotation angle, rotate around Y axis
end_winding_coordinates_rotated(:,1) = end_winding_coordinates_temp(:,1);
end_winding_coordinates_rotated(:,2) = end_winding_coordinates_temp(:,2)*cosd(-end_winding_rotation_angle) + end_winding_coordinates_temp(:,3)*sind(-end_winding_rotation_angle);
end_winding_coordinates_rotated(:,3) = end_winding_coordinates_temp(:,3)*cosd(-end_winding_rotation_angle) - end_winding_coordinates_temp(:,2)*sind(-end_winding_rotation_angle);

%move to z position aligned with the vertical end coils (i.e minimum z of
%the vertical end coil)
end_winding_coordinates_rotated(:,3) = end_winding_coordinates_rotated(:,3) + (min(end_winding_coordinates(:,3))-min(end_winding_coordinates_rotated(:,3)));

%rotate so that aligned with the coils in the stacks
rotation_angle2 = (2*n-1) * machine.pole_angle - machine.pole_angle +0.5 * filament_HTS.coil_angle;

%Rotate initial end winding coordinates
end_winding_coordinates_temp = end_winding_coordinates_rotated;

%rotation angle
end_winding_coordinates_rotated(:,1) = end_winding_coordinates_temp(:,1)*cosd(rotation_angle2) - end_winding_coordinates_temp(:,2)*sind(rotation_angle2);
end_winding_coordinates_rotated(:,2) = end_winding_coordinates_temp(:,2)*cosd(rotation_angle2) + end_winding_coordinates_temp(:,1)*sind(rotation_angle2);
end_winding_coordinates_rotated(:,3) = end_winding_coordinates_temp(:,3);

%Add End Winding Coil
[BSmag] = BSmag_add_filament(BSmag,end_winding_coordinates_rotated,I,dGamma2);

%add the other end coil (mirror over the Z-axis)
end_winding_coordinates_rotated(:,3) = -1 * end_winding_coordinates_rotated(:,3);     % Take symmetry over X-Y plane, by multlipying Z coordinate by -1
[BSmag] = BSmag_add_filament(BSmag,end_winding_coordinates_rotated,I,dGamma2);

end
else
    %do nothing
end
