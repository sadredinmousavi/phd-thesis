% Permanent Magnet Force Field Simulation
clear; close all; clc;

%% Parameters
mu0 = 4*pi*1e-7; % Permeability of free space (T·m/A)
Br = 1.2; % Remanence of the magnet (Tesla) - typical for NdFeB
volume = [0.01, 0.01, 0.01]; % Magnet dimensions [x,y,z] in meters
m = Br*volume(1)*volume(2)*volume(3)/mu0; % Magnetic dipole moment (A·m²)
m_direction = [0, 0, 1]; % Orientation of magnetization (z-direction)

%% Grid Setup
x_range = linspace(-0.05, 0.05, 30); % x-coordinates (m)
y_range = linspace(-0.05, 0.05, 30); % y-coordinates (m)
z_range = linspace(0.01, 0.1, 20);   % z-coordinates (m) (above magnet)

[X, Y, Z] = meshgrid(x_range, y_range, z_range);

%% Calculate Magnetic Field (B) at each point
Bx = zeros(size(X));
By = zeros(size(Y));
Bz = zeros(size(Z));

% Position of the dipole (center of magnet)
r0 = [0, 0, 0];

for i = 1:numel(X)
    % Position vector from dipole to observation point
    r = [X(i), Y(i), Z(i)] - r0;
    distance = norm(r);
    
    if distance > 0
        % Unit vector
        r_hat = r/distance;
        
        % Magnetic field of a dipole (vector form)
        B = (mu0/(4*pi)) * ( (3*dot(m_direction,r_hat)*r_hat - m_direction)/distance^3 );
        
        Bx(i) = B(1);
        By(i) = B(2);
        Bz(i) = B(3);
    end
end

%% Calculate Force on a small test dipole (e.g., another magnet)
test_m = [0, 0, 1e-3]; % Test dipole moment (A·m²)
test_m = test_m/norm(test_m); % Normalize

% Force is gradient of potential energy F = -∇U = ∇(m·B)
[Fx, Fy, Fz] = gradient(-(test_m(1)*Bx + test_m(2)*By + test_m(3)*Bz), ...
                  x_range, y_range, z_range);

%% Visualization
% Select a z-plane to visualize
z_idx = 5; % Index of z-plane to visualize
z_value = z_range(z_idx);

% Magnetic Field
figure(1)
quiver(X(:,:,z_idx), Y(:,:,z_idx), Bx(:,:,z_idx), By(:,:,z_idx))
title(['Magnetic Field at z = ' num2str(z_value) ' m'])
xlabel('x (m)'); ylabel('y (m)');
axis equal; grid on;

% Force Field
figure(2)
quiver(X(:,:,z_idx), Y(:,:,z_idx), Fx(:,:,z_idx), Fy(:,:,z_idx))
title(['Force Field at z = ' num2str(z_value) ' m'])
xlabel('x (m)'); ylabel('y (m)');
axis equal; grid on;

% 3D Visualization
figure(3)
slice(X, Y, Z, sqrt(Bx.^2 + By.^2 + Bz.^2), 0, 0, [0.02 0.05 0.08])
shading interp
colorbar
title('Magnetic Field Magnitude (T)')
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');