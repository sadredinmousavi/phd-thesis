% Simulation of a small magnet moving under the influence of a fixed magnet
clear; close all; clc;

% Physical constants
mu0 = 4*pi*1e-7;   % Permeability of free space (T·m/A)

% Big magnet (fixed at x = -0.1 m)
Br_big = 1.2;      % Remanence (Tesla) - NdFeB magnet
V_big = [0.02, 0.02, 0.02]; % Volume [x,y,z] (m³)
m_big = (Br_big * prod(V_big) / mu0) * [1; 0; 0]; % Dipole moment (A·m², column vector)
pos_big = [-0.1; 0; 0]; % Position (m, column vector)

% Small magnet (movable, starts at x = 0)
Br_small = 1.2;    % Same material
V_small = [0.005, 0.005, 0.005]; % Smaller volume
m_small = (Br_small * prod(V_small) / mu0) * [1; 0; 0]; % Dipole moment (A·m², column vector)
pos_small_initial = [0; 0; 0]; % Initial position (m, column vector)
vel_small_initial = [0; 0; 0]; % Initial velocity (m/s, column vector)

% Simulation parameters
tspan = [0, 5]; % Time span (s)
mass = 0.01;    % Mass of small magnet (kg)

% ODE solver (ode45)
initial_conditions = [pos_small_initial; vel_small_initial]; % Column vector
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, y] = ode45(@(t,y) magnet_dynamics(t, y, m_big, pos_big, m_small, mass, mu0), tspan, initial_conditions, options);

% Extract position and velocity
pos_small = y(:, 1:3);
vel_small = y(:, 4:6);

% Animation
figure;
axis_limit = 0.15;
for k = 1:length(t)
    clf;
    hold on;
    
    % Plot big magnet (red)
    plot3(pos_big(1), pos_big(2), pos_big(3), 'ro', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
    
    % Plot small magnet (blue)
    plot3(pos_small(k,1), pos_small(k,2), pos_small(k,3), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    
    % Plot trajectory (dotted line)
    plot3(pos_small(1:k,1), pos_small(1:k,2), pos_small(1:k,3), 'b:');
    
    % Labels and limits
    xlim([-axis_limit, axis_limit]);
    ylim([-axis_limit, axis_limit]);
    zlim([-axis_limit, axis_limit]);
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
    title(['Time: ', num2str(t(k), '%.2f'), ' s']);
    grid on;
    view(2); % 2D view (x-y plane)
    drawnow;
    
    % Pause for animation speed
    pause(0.05);
end

% Dynamics function for ODE45
function dydt = magnet_dynamics(t, y, m_big, pos_big, m_small, mass, mu0)
    % Extract position and velocity of small magnet (column vectors)
    pos_small = y(1:3);
    vel_small = y(4:6);
    
    % Distance vector between magnets (column vector)
    r = pos_small - pos_big;
    r_norm = norm(r);
    
    % Magnetic field (B) due to big magnet at small magnet's position
    if r_norm > 0
        r_hat = r / r_norm;
        B = (mu0/(4*pi)) * ( (3*(m_big' * r_hat) * r_hat - m_big) / r_norm^3 );
    else
        B = zeros(3,1); % Avoid singularity
    end
    
    % Force on small magnet: F = ∇(m_small · B)
    % Numerical gradient approximation (finite difference)
    eps = 1e-6; % Small perturbation
    F = zeros(3,1);
    for i = 1:3
        r_perturbed = r;
        r_perturbed(i) = r_perturbed(i) + eps;
        r_perturbed_norm = norm(r_perturbed);
        if r_perturbed_norm > 0
            r_perturbed_hat = r_perturbed / r_perturbed_norm;
            B_perturbed = (mu0/(4*pi)) * ( (3*(m_big' * r_perturbed_hat) * r_perturbed_hat - m_big) / r_perturbed_norm^3 );
            F(i) = (m_small' * B_perturbed - m_small' * B) / eps;
        end
    end
    
    % Equations of motion (Newton's 2nd law)
    acceleration = F / mass;
    
    % Return derivatives (column vector)
    dydt = [vel_small; acceleration];
end