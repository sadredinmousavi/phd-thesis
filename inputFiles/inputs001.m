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
eq_times = [10, 20, 30, 40, 50];  % Some times repeat
% eq_times = [10,];
eq_positions = [
     0.05,  0.05;
     0.05,  0.00;  % Duplicate, should merge into one
    -0.05, -0.05;
    -0.03 , 0.03;
     0.0,  0.0
];







eq_times = [0, 5];  % Some times repeat
eq_positions = [
     0.00,  0.00;
     0.02,  0.02;
];

% Generate circular equilibrium points
centerPoint = [0.02, 0.02]; 
r1 = 0.03; 
N = 8; 
s1 = 5;
last_time = eq_times(end);
[eq_positions_1, eq_times_1] = defs.generateCircularEqPoints(centerPoint, r1, N, s1, last_time);


% Generate linear movement points
startPoint = [0.02, 0.02];  
endPoint = [-0.02, -0.02];  
N = 3;  
s1 = 5;
last_time = eq_times_1(end);
[eq_positions_2, eq_times_2] = defs.generateLinearEqPoints(startPoint, endPoint, N, s1, last_time);


% Generate circular equilibrium points
centerPoint = [-0.02, -0.02]; 
r1 = 0.03; 
N = 8; 
s1 = 5;
last_time = eq_times_2(end);
[eq_positions_3, eq_times_3] = defs.generateCircularEqPoints(centerPoint, r1, N, s1, last_time);


last_time = eq_times_3(end);
eq_positions = [eq_positions; eq_positions_1; eq_positions_2; eq_positions_3; 0, 0]; 
eq_times = [eq_times(:); eq_times_1(:); eq_times_2(:); eq_times_3(:); last_time+s1]; 
defs = defs.defineEquilibriumPoints(eq_times, eq_positions);



%%
fieldSim = PermanentMagnetField(defs);
fieldSim.makeForceFunctionHradCode();




%% Run Simulation
% Create a Simulation instance
sim = Simulation(defs);

% Calculate psai_array for all equilibrium points
useSavedData = true; 
filename = 'case_001.mat'; 
sim = sim.calcEquilibrium(useSavedData, filename);
sim.generateExperimentalInputs('experimental_inputs_001.txt');

% snapshot_time = 10;
% idx = find([sim.timed_psai_array.time] == snapshot_time, 1);
% if ~isempty(idx)
%     eq_idx = find([sim.defs.eq_points.time] == snapshot_time, 1);
%     eq_positions = sim.defs.eq_points(eq_idx).positions;
%     sim.plotStaticSnapshot(sim.timed_psai_array(idx).psai_array, eq_positions);
% %     sim.plotStaticSnapshot(sim.timed_psai_array(2).psai_array, sim.defs.eq_points(2).positions); 
% else
%     fprintf('No equilibrium point found for time %.2f\n', snapshot_time);
% end


% for i=1:length(sim.timed_psai_array)
%     sim.plotStaticSnapshot(sim.timed_psai_array(i).psai_array, sim.defs.eq_points(i).positions); 
%     keyboard
% end


tspan = [0 eq_times(end) + 5];
% Initial conditions with zero velocity
initial_positions = sim.defs.microRobots.positions(:,1:2);
initial_velocities = zeros(size(initial_positions));
initial_conditions = [initial_positions, initial_velocities];

% Run simulation with drag
sim = sim.runODESimulation(tspan, initial_conditions);

% Plot results
% sim.plotTrajectories();


% After running your simulation:
sim.createSimulationVideo('robot_simulation_001.mp4', ...
    'ShowForceField', false, ...
    'FrameRate', 3, ...
    'ShowMagnets', true);





