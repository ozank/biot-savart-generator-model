%---------------------------------------------------
%  NAME:      Biot Savart Template.m
%  WHAT:      Template function for the Biot Savart Model that calls the other functions for analyze and plotting
%  AUTHOR:    Ozan Keysan (07/2025)
%  REQUIRED:  BSmag Toolbox 20150407 (for magnetic field calculations)
%----------------------------------------------------

%% Initialize Biot Savart Model
%clear all, close all, clc
BSmag = BSmag_init(); % Initialize BSmag analysis

%BIOT SAVART MODEL SETTINGS
%Discrete steps along coil during, biot savart calculations, has a direct
%effect on calculation time
dGamma2 = 1e-2; % filament max discretization step [m], default to 1 cm  


%% Get Machine Parameters
machine_parameters;      %large machine parameters
%small_machine_parameters;      %large machine parameters

%% Get Material Properties
material_constants; % Load material constants


%% Get Winding coordinates

% It is possible to draw two types of windings: Wave winding and
% conventional race track coils, Please comment out the unwanted type, and
% use ONLY one of the winding types


%% WAVE WINDING
%Get Wave Winding Coordinates
% comment out below for wave winding coil
%wave_winding_coordinates;
%Nmodules_radial = 3; % Number of modules to be simulated in the radial direction, default 3

% Add windings for the wave winding
%plot_wave_winding;

%% RACE TRACK COIL WINDING
%Get Race track axial machine Winding Coordinates 
% comment out below for race track coil

 axial_winding_coordinates;
 Npoles_radial = 4; % Number of modules to be simulated in the radial direction, default 8
 %Npoles_radial = machine.Npole/4; 

% Add windings for the axial race track winding
 plot_axial_race_track_winding;

%% Electrical Machine Parameter Estimations

%Get flux calculations
%calculate_flux_per_pole;  % Outputs flux per pole and maximum B values


%% Biot Savart Calculation Space Definitions
% Field points (where we want to calculate the field)

% % %% Solution Space A: Polar Plane Segment on Z-axis
% % % %Solution Space 
%  R_min = 0.15;  %Solution space inner radius
%  R_max = 0.32; %Solution space outer space
% % 
%  angle_offset = 0; %Solution space starting point
%  angle_span = 360; % Solution angle span (degrees)
% % 
%  data_point_angle= 3600;  % number of data points in the tangential directions (through angle)
%  data_point_radius = 50; %number of data points in the radial (radius) direction
% % 
%  r_M = linspace (R_min,R_max, data_point_radius+1);
%  angle_M = linspace (angle_offset,angle_offset + angle_span, data_point_angle+1);
% % 
% % %create polar coordinate points
%  [R_M,ANGLE_M] = ndgrid(r_M,angle_M * (pi()/180));
% % 
% % %Convert from polar coordinates to cartesian points
%  [X_M,Y_M] = pol2cart(ANGLE_M,R_M);
%  Z_M = zeros(data_point_radius+1,data_point_angle+1); % z [m] 

% Solution SPace B: Cylindrical Surface at a Specific Radius

%Solution Space 
R = HTS.R_mean;  % Radius of the cylindrical surface

Z_min = -0.08;  %Solution space minimum Z point
Z_max = 0.08; %Solution space maximum Z point

angle_offset = 0; %Solution space starting point
angle_span = 45; % Solution angle span (degrees)

data_point_angle= 40;  % number of data points in the tangential directions (through angle)
data_point_radius = 20; %number of data points in the Z direction

r_M = linspace (R, R, data_point_radius+1);  %Constant R values
angle_M = linspace (angle_offset,angle_offset + angle_span, data_point_angle+1);

%create polar coordinate points
[R_M,ANGLE_M] = ndgrid(r_M,angle_M * (pi()/180));

%Convert from polar coordinates to cartesian points
[X_M,Y_M] = pol2cart(ANGLE_M,R_M);

z_M = linspace (Z_min,Z_max, data_point_radius+1);
Z_M = repmat(z_M,  data_point_angle+1,1)';

% BIOT SAVART ANALYSIS RUN
%Biot-Savart Integration
 [BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M);

%% Post Processing 
%BSmag_plot_field_points(BSmag,X_M,Y_M,Z_M); % -> shows the field point line

% Plot B/|B|
figure(1)
    normB=sqrt(BX.^2+BY.^2+BZ.^2);
    %quiver3(X,Y,Z,BX./normB,BY./normB,BZ./normB,'r')
    quiver3(X,Y,Z,BX./normB,BY./normB,BZ./normB,1,'b')
    % hold on
 %contour(X,Z,normB)
 axis equal
% hold off
xlabel ('x [m]'), ylabel ('y [m]'), title ('Bz [T]')
% 

%Plot Bz on the plane
figure(2), hold on, box on, grid on
    contourf(R_M.*sin(ANGLE_M), Z, sqrt(BX.^2+BY.^2+BZ.^2)), colorbar
xlabel ('x [m]'), ylabel ('y [m]'), title ('Bmag [T]')



figure(3), hold on, box on, grid on
    contourf(R_M.*sin(ANGLE_M), Z, sqrt((BX.*sin(ANGLE_M)).^2 + (BY.*cos(ANGLE_M)).^2)), colorbar
    axis equal
xlabel ('x [m]'), ylabel ('y [m]'), title ('B_tangential [T]')


% Plot Bz on the plane
% figure(2), hold on, box on, grid on
%     contourf(X, Y, BZ), colorbar
% xlabel ('x [m]'), ylabel ('y [m]'), title ('Bz [T]')

% %Other Plot Options for Reserve
% 
% % Plot Bz on the volume
% figure(3), hold on, box on, grid on
 Gamma = winding_coordinates_rotated;
% plot3(Gamma(:,1),Gamma(:,2),Gamma(:,3),'.-r') % plot filament
% slice(X,Y,Z,BZ,[0],[],[-1,0,1]), colorbar % plot Bz
% xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Bz [T]')
% view(3), axis equal, axis tight
% caxis([-0.5,0.5]*1e-5)
% 
% 
% Plot some flux tubes

figure(4), hold on, box on, grid on
plot3(Gamma(:,1),Gamma(:,2),Gamma(:,3),'.-r') % plot filament
[X0,Y0,Z0] = ndgrid(-1.5:0.5:1.5,-1.5:0.5:1.5,-2); % define tubes starting point
htubes = streamtube(stream3(X,Y,Z,BX,BY,BZ,X0,Y0,Z0), [0.2 10]);
xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Some flux tubes')
view(3), axis equal, axis tight
set(htubes,'EdgeColor','none','FaceColor','c') % change tube color
camlight left % change tube light

