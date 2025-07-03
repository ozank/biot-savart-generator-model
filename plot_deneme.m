% Field points (where we want to calculate the field)

%vector calculation

%offset Z for different lenghts
x_M = linspace(-0.5,0.5,21); % x [m]
y_M = linspace(-0.5,0.5,21); % y [m]
z_M = linspace(-0.5,0.5,21); % z [m]
%Create mesh grid in polar coordinates

[X_mesh,Y_mesh,Z_mesh]=meshgrid(x_M, y_M , z_M);

BSmag_plot_field_points(BSmag,X_mesh,Y_mesh,Z_mesh); % shows the field points volume
% Biot-Savart Integration
[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_mesh,Y_mesh,Z_mesh);
% Plot B/|B|
figure(1)
 normB=sqrt(BX.^2+BY.^2+BZ.^2);
 quiver3(X,Y,Z,BX./normB,BY./normB,BZ./normB,'b')
%axis tight

% Plot Bz on the volume
figure(2), hold on, box on, grid on
plot3(winding_coordinates_rotated(:,1),winding_coordinates_rotated(:,2),winding_coordinates_rotated(:,3),'.-r') % plot filament
slice(X,Y,Z,BZ,[],[],[-0.1,0,0.1]), colorbar % plot Bz
xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Bz [T]')
view(3), axis equal, axis tight
%caxis([-0.5,0.5]*1e-5)


% Plot some flux tubes
%figure(3), hold on, box on, grid on
%plot3(winding_coordinates_rotated(:,1),winding_coordinates_rotated(:,2),winding_coordinates_rotated(:,3),'.-r') % plot filament
figure;
[X0,Y0,Z0] = ndgrid(x_M,y_M,0); % define tubes starting point
streamlines = stream3(X,Y,Z,BX,BY,BZ,X0,Y0,Z0);
htubes = streamtube(streamlines, [0.25 10]);

xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Some flux tubes')
view(3), axis equal, axis tight
set(htubes,'EdgeColor','none','FaceColor','c') % change tube color
camlight left % change tube light


%stream over a plane at Z coordinate
%streamslice(X,Y,Z,BX,BY,BZ,[],[],0.1)