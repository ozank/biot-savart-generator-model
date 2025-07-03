% Create a new figure
figure;
clf; % Clear the current figure

% Define the 3D grid
[x, y, z] = meshgrid(-2:0.2:2, -2:0.2:2, -2:0.2:2);

% Define a 3D swirling vector field
bx = -y;
by = x;
bz = z;

% Define starting points for streamlines
[x0, y0, z0] = meshgrid(-1:1:1, -1:1:1, -1:1:1);
x0 = x0(:);
y0 = y0(:);
z0 = z0(:);

% Compute streamlines
streamlines = stream3(x, y, z, bx, by, bz, x0, y0, z0);

% Compute vector magnitudes at streamline points (for coloring)
% Concatenate all streamline points
all_points = cell2mat(streamlines');
xi = all_points(:,1);
yi = all_points(:,2);
zi = all_points(:,3);
bi = sqrt(interp3(x, y, z, bx.^2 + by.^2 + bz.^2, xi, yi, zi, 'linear', 0));

% Draw streamtubes with radius 0.2 and 10-sided cross-sections
htubes = streamtube(streamlines, [0.2 10]);

% Get tube vertex color data
if ~isempty(htubes)
    for i = 1:length(htubes)
        set(htubes(i), ...
            'FaceColor', 'interp', ...
            'EdgeColor', 'none', ...
            'FaceLighting', 'gouraud');
    end
end

% Add lighting and shading
camlight('headlight');
lighting gouraud;
shading interp;

% Set axis and view
axis tight
view(3)
grid on
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Streamtubes of a Swirling Vector Field')
