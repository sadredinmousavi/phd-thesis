%% Clear workspace and command window
clear; clc;

%% Create the PermanentMagnetField object and add magnets on a circle.
pmField = PermanentMagnetField;

% Parameters for the magnets
nMagnets = 6;               % number of magnets
r_circle = 0.1;             % circle radius in meters
Br = 1.3;                   % Remanence in Tesla for N42 magnets
dims = [0.02, 0.02, 0.05];    % dimensions of each magnet (m)

% Add n magnets evenly spaced on a circle (XY plane) with radial-in orientation.
pmField = pmField.addMagnetsOnCircle(nMagnets, r_circle, 'radial_in', Br, 'dimensions', dims);

%% Create an observation grid in the XY plane (Try increasing the resolution)
xRange = linspace(-0.25, 0.25, 50);
yRange = linspace(-0.25, 0.25, 50);
[x, y] = meshgrid(xRange, yRange);
z_level = 0;  % observation plane; you might try a nonzero z-level too

% Preallocate arrays for field components.
Bx = zeros(size(x));
By = zeros(size(x));
Bz = zeros(size(x));

% Compute the field at each grid point.
for idx = 1:numel(x)
    [Bx(idx), By(idx), Bz(idx)] = pmField.calculateField(x(idx), y(idx), z_level);
end

% Compute the total field magnitude.
Bmag = sqrt(Bx.^2 + By.^2 + Bz.^2);

%% Visualization

% Figure 1: Quiver plot of field vectors
figure;
% Increase the scale factor (adjust as necessary, try 5 or 10) so that small vectors become visible.
quiver(x, y, Bx, By, 5, 'r', 'LineWidth', 1.5);
xlabel('X (m)');
ylabel('Y (m)');
title('Magnetic Field Vectors in XY Plane (z = 0)');
axis equal; grid on;

% Figure 2: Contour plot of field magnitude
figure;
contourf(x, y, Bmag, 20, 'LineColor', 'none');
colorbar;
xlabel('X (m)');
ylabel('Y (m)');
title('Magnetic Field Magnitude (T) in XY Plane (z = 0)');
axis equal; grid on;


% (Optional) Figure 3: Quiver plot for a different observation plane (e.g., z = 0.02)
z_level_2 = 0.02;
Bx2 = zeros(size(x));
By2 = zeros(size(x));
Bz2 = zeros(size(x));
for idx = 1:numel(x)
    [Bx2(idx), By2(idx), Bz2(idx)] = pmField.calculateField(x(idx), y(idx), z_level_2);
end
figure;
quiver(x, y, Bx2, By2, 5, 'b', 'LineWidth', 1.5);
xlabel('X (m)');
ylabel('Y (m)');
title('Magnetic Field Vectors in XY Plane (z = 0.02)');
axis equal; grid on;
