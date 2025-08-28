% Trial Code for Variable End Winding Angle
% Changes the end winding angle in a foor loop to check its effect

%rotation_angle = 45;

tic
% Move the marker along the parabola and capture frames in a loop
for i=0:10

    HTS.coil_length = 0.1 + 0.05*i;
    biot_savart_template;

  figure(1)
    %normB=sqrt(BX.^2+BY.^2+BZ.^2);
    %quiver3(X,Y,Z,BX./normB,BY./normB,BZ./normB,'r')
      contourf(X, Y, BZ), colorbar
    %quiver3(X,Y,Z,BX./normB,BY./normB,BZ./normB,3,'b')
    % hold on
 %contour(X,Z,normB)
 axis equal
 clim([-3 3])
xlim([1.5,3]); 
 ylim([0,1]); 
zlim([-0.5,0.5]); 
% hold off
xlabel ('x [m]'), ylabel ('y [m]'), title (['Bz airgap [T]',' HTS Coil Length (mm)=', num2str(HTS.coil_length *1000)])
view(-30,40)
exportgraphics(gca,"overview.gif",Append=true)
close;


    % Plot Btangential on the plane
figure(2), hold on, box on, grid on
    contourf(X, Y, BZ), colorbar
    clim([-3 3])
    axis equal
xlim([1.5,3]); 
 ylim([0,1]); 
zlim([-0.5,0.5]); 
xlabel ('x [m]'), ylabel ('y [m]'), title (['Bz airgap [T]',' HTS Coil Length (mm) =', num2str(HTS.coil_length *1000)])
exportgraphics(gca,"B_horizontal.gif",Append=true)
close;

        % Plot Btangential on the plane
% figure(2), hold on, box on, grid on
% contourf(R_M.*sin(ANGLE_M), Z, -BZ), colorbar
% clim([0 5])
% axis equal
% xlabel ('x [m]'), ylabel ('y [m]'), title (['B_Z (Vertical) [T]', 'Coil Angle =', num2str(end_winding_rotation_angle)] )
% exportgraphics(gca,"B_Z.gif",Append=true)
% close;

end

toc

