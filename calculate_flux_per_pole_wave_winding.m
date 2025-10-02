%---------------------------------------------------
%  NAME:      Calculate Flux Per Pole Wave Winding .m
%  WHAT:      For machine sizing problems, calculate the flux per pole for the wave winding
%  The calculations are slightly different as in the wave winding there is an offset in the airgap flux density distribution
% Therefore flux linkage is calculated in two different positions and the induced voltage is calculated according to difference 
% between these two positions

%  AUTHOR:    Ozan Keysan (09/2025)
%----------------------------------------------------
% 
%   Inputs:
%       Main function and motor dimensions
% 
%   Outputs:
% Maximum flux per pole (for wave winding coil) and negative maximum for a second position
% for the positions wave winding is aligned with the stator winding
% The results are used for the induced voltage calculations

%Model to model Single Filament Stator Coil Model
% Stator winding is approximated by single filament passing through the
% mid-point of the stator

%% Stator Filament
%Stator Winding Filament Polar Coordinates(Integral Surface coordinates)
filament.R_inner = stator.R_inner + 0.5 * stator.coil_width;    %Filament modelling Inner Radius
filament.R_outer = stator.R_outer - 0.5 * stator.coil_width;    %Filament modelling Outer Radius

%Stator Coil  Span Angle 
stator.coil_angle = 360 / stator.Ncoil;   % One stator coil pole pitch angle in degrees
stator.coil_pitch = stator.R_mean * (stator.coil_angle *pi() /180); %Coil pitch at mean radius

%Filament length for mean stator coil loop length (for loss calculations)
stator.mean_turn_length = 2*(filament.R_outer - filament.R_inner) + 2 * (stator.coil_pitch - stator.coil_width); %[m] mean turn length of stator coil

%Solution Space Settings
%maximum aligment of the stator coils are shifted by 0.5*machine.pole_angle
%in wave winding type due to geometry difference

angle_offset = 2.5* machine.pole_angle; %Solution space starting point (0 point is the aligned position with the axis)

% Stator dimensions
% Stator Filament Span Angle 
filament.pitch = stator.coil_pitch - stator.coil_width; %Filament pitch at mean radius
filament.span_angle = (filament.pitch / stator.R_mean) * 180/pi(); %Filament pitch angle (degrees)

filament.angle_min = angle_offset + 0.5*(machine.pole_angle - filament.span_angle); %filament starting angle, degrees
filament.angle_max = filament.angle_min + filament.span_angle; %filament finishing angle, degrees

%Careful with the integrations below, if each data point is calculated for
%integration, it slightly overestimates due to extra area covered by dR and
%d_angle calculations below, therefore the solution space should be cropped
%for a more accurate calculation (especially for low number of data points)
%flux density is calculated by Biot savart over data points, 
% integral is calculated using polar coordinated accordingly

%Calculate solution space coordinates
d_R = (filament.R_outer -filament.R_inner)/(data_point_radius -1); % [meters], delta radius, between two data points, required for integral operation
d_ANGLE = filament.span_angle / (data_point_angle -1);  % [degrees], delta angle between two data points, required for integral operation 


r_M = linspace (filament.R_inner + 0.5*d_R ,filament.R_outer-0.5*d_R, data_point_radius-1);      %Radius(m), points for polar coordinates

%First Solution Space
%When the stator coil is aligned with one pole of the wave winding
angle_M1 = linspace (filament.angle_min + 0.5*d_ANGLE, ...
                    filament.angle_max - 0.5*d_ANGLE, ...
                    data_point_angle - 1);   %angle(degrees), points for polar coordinates

%Second solution space
%when the stator coil is aligned with the next wave winding pole
angle_M2 = linspace (machine.pole_angle + filament.angle_min + 0.5*d_ANGLE, ...
                    machine.pole_angle + filament.angle_max - 0.5*d_ANGLE, ...
                    data_point_angle - 1);   %angle(degrees), points for polar coordinates

%Combine these two angle data, to obtain a single solution space
%Biot savart model is run through this combined points
angle_M = [angle_M1, angle_M2];

%create polar coordinate points
[R_M,ANGLE_M] = ndgrid(r_M,angle_M * (pi()/180));

%Convert to rectangular coordinates
[X_M,Y_M] = pol2cart(ANGLE_M,R_M);

%Z coordinate of the filament, can be modified for any axial displacement
Z_M = zeros(size(X_M));%+0.0325; % z [m] %equal size with X_M


BSmag_plot_field_points(BSmag,X_M,Y_M,Z_M); % shows the field points plane

% Biot-Savart Integration
[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M);

%BSmag_plot_field_points(BSmag,X_M,Y_M,Z_M); % -> shows the field point line

%% Calculate Flux Linkage per coil
% (Bz, surface integral over polar coordinates)
%https://math.stackexchange.com/questions/145939/simple-proof-of-integration-in-polar-coordinates

% As there is an offset in the airgap flux density,
%flux linkage should be calculated in two positions where the stator coils
%are aligned with the winding. The difference between these two positions/2
%is defined as the flux per pole, which is used for induced voltage
%calculations

% The BZ data is divided into two halves and flux per pole is calculated
% for each half separately. 
data_point_half = size(BZ,2)/2; %find the data size

% Positive cycle of the stator - pole alignment
stator.flux_per_pole_positive = sum (BZ(:,1:data_point_half) .* R_M(:,1:data_point_half) *d_R * d_ANGLE*(pi/180), "all");   %[Wb], Maximum flux in the stator coils
% int (Bz.dA) in polar coordinates,

% Negative cycle of the stator - pole alignment
stator.flux_per_pole_negative = sum (BZ(:,data_point_half+1 : end) .* R_M(:,data_point_half+1 : end) *d_R * d_ANGLE*(pi/180), "all");   %[Wb], Maximum flux in the stator coils
% int (Bz.dA) in polar coordinates,

%For induced voltage calculations
% flux per pole is calculated as ((flux_positive)-(flux_negative))/2
% this is to reduce the offset and find the flux change in one pole pair
% movement
stator.flux_per_pole = ( stator.flux_per_pole_positive - stator.flux_per_pole_negative)/2; %equivalent flux per pole


%% Maximum Airgap Flux Density
%Maximum flux density in the airgap to calculate the eddy current losses in
%the stator, Gethe the average of the top 10% values to prevent issues with
%calculation singularities 

%positive flux density values
%Get the half of the data set
temp = BZ(:,1:data_point_half);

stator.B_max_positive = mean(maxk(temp(:), ceil(numel(temp)*0.1))); %[Tesla], maximum airgap flux density, averaged over top %10 values

%negative flux density values
temp = BZ(:, data_point_half+1 : end);

stator.B_max_negative =  mean(mink(temp(:), ceil(numel(temp)*0.1))); %[Tesla], maximum airgap flux density, averaged over top %10 values

%Similar to flux per pole calculations for AC loss calculations
%The B_max is calculated using the difference between positive and negative
%cycles (divided by 2 for averaging)

stator.B_max = (stator.B_max_positive - stator.B_max_negative)/2;
%stator.B_max = mean(maxk(BZ(:), ceil(numel(BZ)*0.1))); %[Tesla], maximum airgap flux density, averaged over top %10 values

% Plot Bz on the plane

%Only plot figures if plot_figures = 1 
%if the variable is not defined it the code will generate figures by
%default plot_figures can be set to 0 to reduce computation time

global plot_figures

if   exist('plot_figures','var') && isscalar(plot_figures) && ~plot_figures
    %Do not plot figures
  
else

    %Plot Figures

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

    %to draw boundaries of the stator windings
%    plot(polyshape([filament.R_inner*cosd(filament.angle_min) filament.R_inner*cosd(filament.angle_max) filament.R_outer*cosd(filament.angle_max) filament.R_outer*cosd(filament.angle_min) ] ...
%    ,[filament.R_inner*sind(filament.angle_min) filament.R_inner*sind(filament.angle_max) filament.R_outer*sind(filament.angle_max) filament.R_outer*sind(filament.angle_min)  ]));

end

