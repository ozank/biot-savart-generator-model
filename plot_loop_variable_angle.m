% Trial Code for Variable End Winding Angle
% Changes the end winding angle in a foor loop to check its effect

%rotation_angle = 45;

tic
% Move the marker along the parabola and capture frames in a loop
for i=0:10

    end_winding_rotation_angle=85+i;
    main;

  figure(1)
    normB=sqrt(BX.^2+BY.^2+BZ.^2);
    %quiver3(X,Y,Z,BX./normB,BY./normB,BZ./normB,'r')
    quiver3(X,Y,Z,BX./normB,BY./normB,BZ./normB,3,'b')
    % hold on
 %contour(X,Z,normB)
 axis equal
% hold off
xlabel ('x [m]'), ylabel ('y [m]'), title (['B vector [T]','Coil Angle =', num2str(end_winding_rotation_angle)])
view(-90,5)
exportgraphics(gca,"overview.gif",Append=true)
close;


    % Plot Btangential on the plane
figure(2), hold on, box on, grid on
    contourf(R_M.*sin(ANGLE_M), Z, sqrt((BX.*sin(ANGLE_M)).^2 + (BY.*cos(ANGLE_M)).^2)), colorbar
    clim([0 3])
    axis equal
xlabel ('x [m]'), ylabel ('y [m]'), title (['|B_{tangential}| (Horizontal) [T]', 'Coil Angle =', num2str(end_winding_rotation_angle)] )
exportgraphics(gca,"B_horizontal.gif",Append=true)
close;

        % Plot Btangential on the plane
figure(2), hold on, box on, grid on
contourf(R_M.*sin(ANGLE_M), Z, -BZ), colorbar
clim([0 5])
axis equal
xlabel ('x [m]'), ylabel ('y [m]'), title (['B_Z (Vertical) [T]', 'Coil Angle =', num2str(end_winding_rotation_angle)] )
exportgraphics(gca,"B_Z.gif",Append=true)
close;

end

toc

