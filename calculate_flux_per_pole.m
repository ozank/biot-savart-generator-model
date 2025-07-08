%Calculate Flux Per Pole

stator.R_outer = HTS.R_outer;  %to be changed
stator.R_inner = HTS.R_inner;  %to be changed
stator.R_mean = HTS.R_mean;  %to be changed


%Model to model Single Filament Stator Coil Model

%Stator Winding Filament Polar Coordinates(Integral Surface coordinates)
filament.R_inner = stator.R_inner + 0.5 * stator.coil_width;    %Filament modelling Inner Radius
filament.R_outer = stator.R_outer - 0.5 * stator.coil_width;    %Filament modelling Outer Radius

%Stator Coil  Span Angle 
stator.coil_angle = 360 / machine.Ncoil;   % One stator coil pole pitch angle in degrees
stator.coil_pitch = stator.R_mean * (stator.coil_angle *pi() /180); %Coil pitch at mean radius


%Stator Stator Filament Span Angle 
filament.pitch = stator.coil_pitch - stator.coil_width; %Filament pitch at mean radius
filament.span_angle = (filament.pitch / stator.R_mean) * 180/pi(); %Filament pitch angle (degrees)

%Solution Space Settings 
angle_offset = machine.pole_angle; %Solution space starting point (0 point is the aligned position with the first pole)

data_point_angle= 10;  % number of data points in the tangential directions (through angle)
data_point_radius = 20; %number of data points in the radial (radius) direction

%Calculate solution space coordinates
r_M = linspace (filament.R_inner,filament.R_outer, data_point_radius);      %Radius range in polar coordinates

angle_M = linspace (angle_offset + 0.5*(machine.pole_angle - filament.span_angle), ...
                    angle_offset + 0.5*(machine.pole_angle + filament.span_angle), ...
                    data_point_angle);   %angle range in polar coordinates

%create polar coordinate points
[R_M,ANGLE_M] = ndgrid(r_M,angle_M * (pi()/180));

%Convert to rectangular coordinates
[X_M,Y_M] = pol2cart(ANGLE_M,R_M);

%Z coordinate of the filament, can be modified for any axial displacement
Z_M = zeros(data_point_radius,data_point_angle); % z [m] %data sirasini kontrol et



%BSmag_plot_field_points(BSmag,X_M,Y_M,Z_M); % shows the field points plane
tic
% Biot-Savart Integration
[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M);
toc
%BSmag_plot_field_points(BSmag,X_M,Y_M,Z_M); % -> shows the field point line


