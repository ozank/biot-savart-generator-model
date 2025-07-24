%---------------------------------------------------
%  NAME:      Calculate Flux Per Pole .m
%  WHAT:      For machine sizing problems, calculate the flux per pole according to 
% machine parameters, will be used for induced voltage calculations for the stator
%  AUTHOR:    Ozan Keysan (07/2025)
%----------------------------------------------------
% 
%   Inputs:
%       Main function and motor dimensions
% 
%   Outputs:
%       Maximum flux per pole for one stator coil (when stator coil is aligned
% with the HTS field winding
%

% Call the main function to draw the coils according to machine dimensions
main;

%% Modify Stator Dimensions (can be kept same)
stator.R_outer = HTS.R_outer+0.1;  %to be changed
stator.R_inner = HTS.R_inner-0.1;  %to be changed
stator.R_mean = HTS.R_mean;  %to be changed



%Model to model Single Filament Stator Coil Model
% Stator winding is approximated by single filament passing through the
% mid-point of the stator

%Stator Winding Filament Polar Coordinates(Integral Surface coordinates)
filament.R_inner = stator.R_inner + 0.5 * stator.coil_width;    %Filament modelling Inner Radius
filament.R_outer = stator.R_outer - 0.5 * stator.coil_width;    %Filament modelling Outer Radius

%Stator Coil  Span Angle 
stator.coil_angle = 360 / stator.Ncoil;   % One stator coil pole pitch angle in degrees
stator.coil_pitch = stator.R_mean * (stator.coil_angle *pi() /180); %Coil pitch at mean radius

%Filament length for mean stator coil lenght (for loss calculations)
stator.mean_turn_length = 2*(filament.R_outer - filament.R_inner) + 2 * (stator.coil_pitch - stator.coil_width); %[m] maen turn length of stator coil

%Solution Space Settings 
angle_offset = 2* machine.pole_angle; %Solution space starting point (0 point is the aligned position with the axis)

data_point_angle= 20;  % number of data points in the tangential directions (through angle)
data_point_radius = 50; %number of data points in the radial (radius) direction

%Stator Stator Filament Span Angle 
filament.pitch = stator.coil_pitch - stator.coil_width; %Filament pitch at mean radius
filament.span_angle = (filament.pitch / stator.R_mean) * 180/pi(); %Filament pitch angle (degrees)

filament.angle_min = angle_offset + 0.5*(machine.pole_angle - filament.span_angle); %filament starting angle, degrees
filament.angle_max = filament.angle_min + filament.span_angle; %filament finishing angle, degrees

%Careful with the integrations below, if each data point is calculated for
%integration, it slightly overestimates due to extra area covered by dR and
%d_angle calculations below, therefore the solution space should be cropped
%for a more accurate calculation (especially for low number of data points)
%Biot savart is made over data points, integral should be calculated
%accordingly

%Calculate solution space coordinates
d_R = (filament.R_outer -filament.R_inner)/(data_point_radius -1); % [meters], delta radius, between two data points, required for integral operation
d_ANGLE = filament.span_angle / (data_point_angle -1);  % [degrees], delta angle between two data points, required for integral operation 


r_M = linspace (filament.R_inner + 0.5*d_R ,filament.R_outer-0.5*d_R, data_point_radius-1);      %Radius(m), points for polar coordinates

%cropped solution space
angle_M = linspace (filament.angle_min + 0.5*d_ANGLE, ...
                    filament.angle_max - 0.5*d_ANGLE, ...
                    data_point_angle-1);   %angle(degrees), points for polar coordinates

%create polar coordinate points
[R_M,ANGLE_M] = ndgrid(r_M,angle_M * (pi()/180));

%Convert to rectangular coordinates
[X_M,Y_M] = pol2cart(ANGLE_M,R_M);

%Z coordinate of the filament, can be modified for any axial displacement
Z_M = zeros(size(X_M)); % z [m] %equal size with X_M


BSmag_plot_field_points(BSmag,X_M,Y_M,Z_M); % shows the field points plane

% Biot-Savart Integration
[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M);

%BSmag_plot_field_points(BSmag,X_M,Y_M,Z_M); % -> shows the field point line

%% Calculate Flux Linkage
% (Bz, surface integral over polar coordinates)
%https://math.stackexchange.com/questions/145939/simple-proof-of-integration-in-polar-coordinates

flux = sum (BZ .* R_M *d_R * d_ANGLE*(pi/180), "all")   % int (Bz.dA) in polar coordinates

% Plot Bz on the plane
figure(2), hold on, box on, grid on
contourf(X,Y, BZ), colorbar
%clim([0 5])
axis equal
xlabel ('x [m]'), ylabel ('y [m]'), title ('Bz [T]')
%add data points
hold on
plot(X(:),Y(:),'k.')
%plot original filament (middle of stator coils area)
hold on
plot(polyshape([filament.R_inner*cosd(filament.angle_min) filament.R_inner*cosd(filament.angle_max) filament.R_outer*cosd(filament.angle_max) filament.R_outer*cosd(filament.angle_min) ] ...
    ,[filament.R_inner*sind(filament.angle_min) filament.R_inner*sind(filament.angle_max) filament.R_outer*sind(filament.angle_max) filament.R_outer*sind(filament.angle_min)  ]));

