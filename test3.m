clear; close all; clc;

% Physical constants
mu0 = 4 * pi * 1e-7; % Permeability of free space (T·m/A)

% Magnet parameters
Br = 1.2; % Remanence (Tesla) - NdFeB magnet
V_magnet = [0.01, 0.01, 0.01]; % Volume [x,y,z] (m³)
m_mag = Br * prod(V_magnet) / mu0; % Magnetic moment magnitude (A·m²)

% Positions of four magnets (equally spaced around a circle, radius = 0.1 m)
angles = linspace(0, 2*pi, 5); % 5 angles (0°, 90°, 180°, 270°)
angles = angles(1:4); % Take first 4 (avoid duplicate at 360°)
radius = 0.1; % Distance from center (m)

% Initialize magnet positions and dipole moments (tilted 25° toward center)
magnet_positions = zeros(4, 3); % [x, y, z]
magnet_moments = zeros(4, 3);   % [mx, my, mz]

for k = 1:4
    % Position of k-th magnet
    x = radius * cos(angles(k));
    y = radius * sin(angles(k));
    magnet_positions(k, :) = [x, y, 0];
    
    % Dipole moment direction (tilted 25° toward center)
    tilt_angle = 55 * pi / 180; % 25° in radians
    m_dir = -[x, y, 0] / radius; % Points toward center (normalized)
    m_dir_tilted = [m_dir(1), m_dir(2), sin(tilt_angle)]; % Tilted 25° upward
    m_dir_tilted = m_dir_tilted / norm(m_dir_tilted); % Normalize
    
    magnet_moments(k, :) = m_mag * m_dir_tilted;
end

% Grid setup (2D plane, z=0)
x_range = linspace(-0.15, 0.15, 30);
y_range = linspace(-0.15, 0.15, 30);
[X, Y] = meshgrid(x_range, y_range);
Z = zeros(size(X));

% Compute B-field and force on a test dipole (aligned with z-axis)
test_m = [0, 0, 1e-3]; % Test dipole moment (A·m²)
Bx = zeros(size(X));
By = zeros(size(X));
Bz = zeros(size(X));

for i = 1:numel(X)
    r_test = [X(i), Y(i), Z(i)]; % Test position
    
    % Total B-field from all four magnets
    B_total = [0, 0, 0];
    for k = 1:4
        r = r_test - magnet_positions(k, :);
        r_norm = norm(r);
        
        if r_norm > 0
            r_hat = r / r_norm;
            m = magnet_moments(k, :);
            
            % Dipole field formula
            B = (mu0 / (4 * pi)) * ( (3 * dot(m, r_hat) * r_hat - m) / r_norm^3 );
            B_total = B_total + B;
        end
    end
    
    Bx(i) = B_total(1);
    By(i) = B_total(2);
    Bz(i) = B_total(3);
end

% Force on test dipole: F = ∇(m_test · B)
[Fx, Fy] = gradient(-test_m(3) * Bz, x_range, y_range);

% Normalize force for better visualization
F_magnitude = sqrt(Fx.^2 + Fy.^2);
Fx_normalized = Fx ./ (F_magnitude + 1e-10); % Avoid division by zero
Fy_normalized = Fy ./ (F_magnitude + 1e-10);

% Plotting
figure;
hold on;

% Plot magnet positions (red circles)
for k = 1:4
    plot(magnet_positions(k, 1), magnet_positions(k, 2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
end

% Quiver plot of force field (normalized)
quiver(X, Y, Fx_normalized, Fy_normalized, 0.5, 'b'); % 0.5 = arrow scaling

% Labels and title
xlabel('x (m)');
ylabel('y (m)');
title('2D Force Field from Four Magnets (Tilted 25°)');
axis equal;
grid on;

% Add a center point
plot(0, 0, 'kx', 'MarkerSize', 12, 'LineWidth', 2);