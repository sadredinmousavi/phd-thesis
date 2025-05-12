%% Defining parameters
 % all paramers have to be set in here
 % s 

clc, close all, clear all
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
cd('..');
addpath(genpath('class'));


%% Define parameters
% Create an instance of the Definitions class.
defs = Definitions();

% Set parameters: 8 magnets arranged on a circle of radius 0.25 m,
numMagnets = 8;
radius = 0.25;
defs = defs.defineExternalMagnets(numMagnets, radius, 0.02, 0.02, 0.02, 'mu_0', 4*pi*1e-7, 'M', 868e3, 'Br', 1.3);

% Define observation space: XY plane, z = 0, resolution = 50 points per axis
defs = defs.setObservationSpace([-0.3, 0.3], [-0.3, 0.3], 0, 50);

% Define groups of micro robots
groupParams = [
    struct('center', [0.01, 0.01], 'gridSize', [3, 3], 'spacing', [0.05, 0.05]), ...
%     struct('center', [-0.1, -0.1], 'gridSize', [2, 4], 'spacing', [0.04, 0.04])
];
defaultMagnetParams = struct( ...
    'dimensions', [0.0005, 0.0005, 0.0005], ... % [a, b, c] in meters
    'M', 1000e3, ... % Magnetization [A/m] n35
    'rho', 7.5e3, ... % Density of material [kg/m^3]
    'Mu', 0.001, ... % Dynamic viscosity [Pa·s] (default for water)
    'orientation', [0, 0, 1] ...
);

% Add micro robots to Definitions
defs = defs.defineMicroRobots(groupParams, defaultMagnetParams);

% Define equilibrium points at times t = 1, 2, and 3 seconds
eq_times = [10, 10, 20, 20, 30];  % Some times repeat
eq_positions = [
     0.05,  0.05;
     0.05,  0.05;  % Duplicate, should merge into one
    -0.02, -0.02;
     0.03 , 0.03;
     0.0,  0.0
];

defs = defs.defineEquilibriumPoints(eq_times, eq_positions);

% Check stored equilibrium points
disp(defs.eq_points);


%%
fieldSim = PermanentMagnetField(defs);
fieldSim.makeForceFunctionHradCode




%% Run Simulation
% Create a Simulation instance
sim = Simulation(defs);


% Calculate psai_array for all equilibrium points
sim = sim.calcEquilibrium();

% Run the simulation for a total duration of 5 seconds
% sim.runSimulation(5);


snapshot_time = 10;  % The time you want to visualize
idx = find([sim.timed_psai_array.time] == snapshot_time, 1);

if ~isempty(idx)
    eq_idx = find([sim.defs.eq_points.time] == snapshot_time, 1);
    
    eq_positions = sim.defs.eq_points(eq_idx).positions;
    
     sim.plotStaticSnapshot(sim.timed_psai_array(idx).psai_array, eq_positions);
else
    fprintf('No equilibrium point found for time %.2f\n', snapshot_time);
end






% Assuming you already have a Simulation object 'sim'
tspan = [0 30]; % 10 second simulation
% Initial conditions with zero velocity
initial_positions = sim.defs.microRobots.positions(:,1:2);
initial_velocities = zeros(size(initial_positions));
initial_conditions = [initial_positions, initial_velocities];

% Run simulation with drag
sim = sim.runODESimulation(tspan, initial_conditions);

% Plot results
% sim.plotTrajectories();


% After running your simulation:
sim.createSimulationVideo('robot_simulation.mp4', ...
    'FrameRate', 30);




%%

% Retrieve a microrobot's properties (from your Definitions.microRobots.flat list):
P_robot = defs.microRobots.positions(1, :);  % e.g., first microrobot's position
% P_robot = [0,0, 0];
m_robot = defs.microRobots.moments(1, :);         % corresponding dipole moment

% Calculate the net force on this microrobot:
F_net = fieldSim.calculateForceOnMicrorobot(P_robot, m_robot);



%% Force Field Visualization in Observation Area

% Retrieve observation grid from defs
X = defs.observation.grid.X;
Y = defs.observation.grid.Y;
Z = defs.observation.grid.Z;

% Initialize force arrays for each grid point
Fx = zeros(size(X));
Fy = zeros(size(Y));
Fz = zeros(size(Z));

% Loop over each grid point and calculate force
for xIdx = 1:size(X, 1)
    for yIdx = 1:size(X, 2)
        % Microrobot's position at this grid point
        P_robot = [X(xIdx, yIdx), Y(xIdx, yIdx), Z(xIdx, yIdx)];
        
        % Assume microrobot has a default magnetic moment
        m_robot = defaultMagnetParams.M * prod(defaultMagnetParams.dimensions) * defaultMagnetParams.orientation;
        
        % Calculate the force on the microrobot at this position
        F = fieldSim.calculateForceOnMicrorobotInPlane(P_robot, m_robot);
        
        % Store the force components in the arrays
        Fx(xIdx, yIdx) = F(1);
        Fy(xIdx, yIdx) = F(2);
        Fz(xIdx, yIdx) = F(3); % Optional: Use if you need Z-component visualization
    end
end

%% Visualize the Force Field
figure;
F_magnitude = sqrt(Fx.^2 + Fy.^2); % Magnitudes of force vectors
% Example for different threshold ranges
thresholds = [1e3, 1e5; 1e5, 1e7; 1e7, inf]; % Rows define [low, high]
for t = 1:size(thresholds, 1)
    % Apply specific range
    mask = (F_magnitude >= thresholds(t, 1)) & (F_magnitude <= thresholds(t, 2));
    Fx_filtered = Fx .* mask;
    Fy_filtered = Fy .* mask;

    % Normalize and plot
    Fx_normalized = Fx_filtered ./ F_magnitude;
    Fy_normalized = Fy_filtered ./ F_magnitude;
    Fx_normalized(~mask) = 0;
    Fy_normalized(~mask) = 0;

%     % Create quiver plot
%     figure;
%     quiver(X, Y, Fx_normalized, Fy_normalized, 'AutoScale', 'on', 'LineWidth', 1.5);
%     title(sprintf('Force Field with Threshold [%d, %d]', thresholds(t, 1), thresholds(t, 2)));
%     xlabel('X');
%     ylabel('Y');
%     axis equal;
%     grid on;
end





%%
%% Compare Force Calculations Across the Workspace

% Retrieve observation grid from defs
X = defs.observation.grid.X;
Y = defs.observation.grid.Y;
Z = defs.observation.grid.Z; % Should be constant at z = 0

% Initialize arrays to store differences
diff_Fx = zeros(size(X));
diff_Fy = zeros(size(Y));
diff_Fz = zeros(size(Z)); % Optional, should ideally be zero since Fz = 0 in-plane

% Loop over the entire workspace grid
for xIdx = 1:size(X, 1)
    for yIdx = 1:size(X, 2)
        % Microrobot's position at this grid point
        P_robot = [X(xIdx, yIdx), Y(xIdx, yIdx), Z(xIdx, yIdx)];

        % Assume microrobot has a default magnetic moment
        m_robot = defaultMagnetParams.M * prod(defaultMagnetParams.dimensions) * defaultMagnetParams.orientation;

        % Calculate force using the original function
        F_orig = fieldSim.calculateForceOnMicrorobot(P_robot, m_robot);

        % Calculate force using the hard-coded function
        x = P_robot(1);
        y = P_robot(2);
        psai_array = defs.epms.psai; % Magnet rotation angles
        F_fast = force_field_fast(x, y, psai_array);

        % Compute and store the difference in forces
        diff_Fx(xIdx, yIdx) = abs(F_orig(1) - F_fast(1));
        diff_Fy(xIdx, yIdx) = abs(F_orig(2) - F_fast(2));
        diff_Fz(xIdx, yIdx) = abs(F_orig(3) - F_fast(3)); % Optional
    end
end

%% Analyze and Visualize Differences

% Compute the norm of the differences for visualization
diff_norm = sqrt(diff_Fx.^2 + diff_Fy.^2 + diff_Fz.^2);

% Plot the differences
figure;
surf(X, Y, log10(diff_norm), 'EdgeColor', 'none');
colorbar;
title('Log10 of Differences Between Original and Hard-coded Force Calculations');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Log10(Difference Norm)');
view(2); % Top-down view

