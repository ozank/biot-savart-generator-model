%---------------------------------------------------
%  NAME:      Plot Axial Race Track Winding.m
%  WHAT:      Plots the conventional race track winding and calculates the magnetic field using the Biot-Savart Calculation Method
%  AUTHOR:    Ozan Keysan (06/2025)
%  REQUIRED:  BSmag Toolbox 20150407 (for magnetic field calculations)
%----------------------------------------------------

% Initialize Biot Savart Model
clear all, close all, clc
BSmag = BSmag_init(); % Initialize BSmag analysis

%Get Machine Parameters
%machine_parameters;
small_machine_parameters;

%Get Wave Winding Coordinates
axial_winding_coordinates;



%BIOT SAVART MODEL SETTINGS

%Define the total current for biot savart conductor
I = HTS.current * HTS.N_turns * HTS.N_layers; %[A] filament current

%Discrete steps along coil during, biot savart calculations
dGamma2 = 1e-2; % filament max discretization step [m]  



Npoles_radial =32; % Number of modules to be simulated in the radial direction, default 8

%Create stack coils 
% !!!
%needs to be defined according to independent machine parameters!
coil_to_coil_gap = 0.1; %ignoring coil thicknes, mid plane of the coils should be considered

% coils are placed symmetrically around Z-axis
Coil_z_min = -0.5 *coil_to_coil_gap * machine.Nstacks %minimum coil z(-) for the outer stack



for s = 1:machine.Nstacks+1 %Create machine.Nstacks+1 number of coils in the axial direction, midplane Z=0

winding_coordinates_stack = winding_coordinates;

winding_coordinates_stack(:,3) = Coil_z_min + (s-1)*coil_to_coil_gap; %z coordinate of -z stack

%[BSmag] = BSmag_add_filament(BSmag,winding_coordinates_stack,I2,dGamma2);


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


% Field points (where we want to calculate the field)

%Solution Space 
R_min = 0;  %Solution space inner radius
R_max = 1; %Solution space outer space

angle_offset = 0; %Solution space starting point
angle_span = 360; % Solution angle span (degrees)

data_point_angle= 180;  % number of data points in the tangential directions (through angle)
data_point_radius = 20; %number of data points in the radial (radius) direction

r_M = linspace (R_min,R_max, data_point_radius+1);
angle_M = linspace (angle_offset,angle_offset + angle_span, data_point_angle+1);

%create polar coordinate points
%[X_M,Y_M] = ndgrid(x_M,y_M);
[R_M,ANGLE_M] = ndgrid(r_M,angle_M * (pi()/180));


[X_M,Y_M] = pol2cart(ANGLE_M,R_M);
Z_M = zeros(data_point_radius+1,data_point_angle+1); % z [m] %data sirasini kontrol et



%BSmag_plot_field_points(BSmag,X_M,Y_M,Z_M); % shows the field points plane
tic
% Biot-Savart Integration
[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M);
toc
BSmag_plot_field_points(BSmag,X_M,Y_M,Z_M); % -> shows the field point line


% Plot B/|B|
figure(1)
    normB=sqrt(BX.^2+BY.^2+BZ.^2);
    %quiver3(X,Y,Z,BX./normB,BY./normB,BZ./normB,'r')
    quiver3(X,Y,Z,BX,BY,BZ,'r')

% Plot Bz on the plane
figure(2), hold on, box on, grid on
    contourf(X, Y, BZ), colorbar
xlabel ('x [m]'), ylabel ('y [m]'), title ('Bz [T]')


% Field points (where we want to calculate the field)
x_M = linspace(-0.5,0.5,21); % x [m]
y_M = linspace(-0.5,0.5,22); % y [m]
z_M = linspace(-0.5,0.5,23); % z [m]
[X_M,Y_M,Z_M]=meshgrid(x_M,y_M,z_M);
BSmag_plot_field_points(BSmag,X_M,Y_M,Z_M); % shows the field points volume
% Biot-Savart Integration
[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M);

% Plot Bz on the volume
figure(3), hold on, box on, grid on
Gamma = winding_coordinates_rotated;
plot3(Gamma(:,1),Gamma(:,2),Gamma(:,3),'.-r') % plot filament
slice(X,Y,Z,BZ,[0],[],[-1,0,1]), colorbar % plot Bz
xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Bz [T]')
view(3), axis equal, axis tight
caxis([-0.5,0.5]*1e-5)


% Plot some flux tubes
%koordinatlar girilecek
figure(4), hold on, box on, grid on
plot3(Gamma(:,1),Gamma(:,2),Gamma(:,3),'.-r') % plot filament
[X0,Y0,Z0] = ndgrid(-1.5:0.5:1.5,-1.5:0.5:1.5,-2); % define tubes starting point
htubes = streamtube(stream3(X,Y,Z,BX,BY,BZ,X0,Y0,Z0), [0.2 10]);
xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Some flux tubes')
view(3), axis equal, axis tight
set(htubes,'EdgeColor','none','FaceColor','c') % change tube color
camlight left % change tube light



