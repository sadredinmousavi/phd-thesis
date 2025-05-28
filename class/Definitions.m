classdef Definitions
    properties
        epms % Struct to store magnet data
        observation   % Stores observation space parameters
        microRobots   % Stores micro robot groups and their properties
        eq_points % Struct to store equilibrium points and associated times
    end
    
    methods
        % Constructor to initialize epms struct
        function obj = Definitions()
            obj.epms = struct('inputs', [], 'MagPos', [], 'psai', [], 'm_vector', []);
            obj.observation = struct('x_space', [], 'y_space', [], 'z_obs', [], 'grid', []);
            obj.microRobots = struct('positions', [], 'moments', [], 'masses', []);
            obj.eq_points = struct('time', [], 'position', []); % Store equilibrium point times and [x, y] positions
        end
        
        % Method to define external magnets
        function obj = defineExternalMagnets(obj, numMagnets, radius, magnetLength, magnetWidth, magnetHeight, varargin)
            % Inputs:
            %   numMagnets     - Number of magnets
            %   radius         - Circular arrangement radius [m]
            %   magnetLength   - Length of the rectangular magnet [m]
            %   magnetWidth    - Width of the rectangular magnet [m]
            %   magnetHeight   - Height of the rectangular magnet [m]
            % Optional inputs (name-value pairs):
            %   'mu_0'         - Override permeability
            %   'M'            - Override magnetization
            
            % Default physical constants
            mu_0 = 4 * pi * 1e-7; % Vacuum permeability [H/m]
            Br = 1.3;             % Magnet Remanence Level [T]
            M = 868e3;            % Magnetization [A/m]
            Mu = 0.001;           % Dynamic viscosity [Pa·s]
            
            % Parse optional inputs
            for i = 1:2:length(varargin)
                switch varargin{i}
                    case 'mu_0'
                        mu_0 = varargin{i+1}; % Override mu_0
                    case 'M'
                        M = varargin{i+1}; % Override M
                    case 'Br'
                        Br = varargin{i+1}; % Override Br
                    case 'Mu'
                        Mu = varargin{i+1}; % Override Mu
                end
            end
            
            % Store inputs in epms.inputs
            obj.epms.inputs = struct( ...
                'numMagnets', numMagnets, ...
                'radius', radius, ...
                'magnetLength', magnetLength, ...
                'magnetWidth', magnetWidth, ...
                'magnetHeight', magnetHeight, ...
                'mu_0', mu_0, ...
                'Br', Br, ...
                'M', M ...
            );
            
            % Calculate volume and magnetic moment
            V = magnetLength * magnetWidth * magnetHeight; % Volume [m^3]
            m = M * V; % Magnetic moment [A*m^2]
            obj.epms.inputs.V = V;
            obj.epms.inputs.m = m;
            
            % Define angular positions (phi) and rotation angles (psai)
            phi = linspace(0, 2*pi, numMagnets + 1); % Evenly spaced angles
            phi = phi(1:end-1);                     % Remove overlap at 2*pi
            psai = zeros(1, numMagnets);            % Initialize rotation angles
            
            % Calculate positions and store them in epms.MagPos
            MagPos = zeros(numMagnets, 3); % [x, y, z]
            for i = 1:numMagnets
                MagPos(i, :) = [radius * cos(phi(i)), radius * sin(phi(i)), 0]; % Position in 2D workspace
            end
            
            % Store in epms
            obj.epms.positions = MagPos;
            obj.epms.psai = psai;
            
            % Calculate magnetic moment vectors and store them in epms.m_vector
            obj.epms.orientations = zeros(numMagnets, 3); % Initialize magnetic moment vectors
            for i = 1:numMagnets
                psai_ = psai(i); % Rotation angle for each magnet
                phi_ = phi(i);
                obj.epms.orientations(i, :) = obj.calculateMagneticMoment(psai_, phi_,  1); % Calculate magnetic moment
            end
        end
        
        % Method to calculate magnetic moment vector dynamically
        function m_vector = calculateMagneticMoment(~, psi, phi, m)
            % Compute [m_x, m_y, m_z] based on psi (rotation about axis)
            m_vector = m * [-sin(phi) * sin(psi), cos(phi) * sin(psi), cos(psi)];
        end
        
        
        
        % Method to define observation-related parameters
        function obj = setObservationSpace(obj, xRange, yRange, zLevel, resolution)
            % Define observation space parameters
            % Inputs:
            %   xRange     - Range for x-axis [xmin, xmax]
            %   yRange     - Range for y-axis [ymin, ymax]
            %   zLevel     - Fixed z-level for observation plane
            %   resolution - Number of points along each axis
            
            x_space = linspace(xRange(1), xRange(2), resolution); % Generate x-axis grid
            y_space = linspace(yRange(1), yRange(2), resolution); % Generate y-axis grid
            z_obs = zLevel; % Fixed z-level for the observation plane
            
            [X, Y] = meshgrid(x_space, y_space); % Create 2D observation grid
            Z = z_obs * ones(size(X)); % Observation plane at fixed z
            
            % Store observation space parameters
            obj.observation.x_space = x_space;
            obj.observation.y_space = y_space;
            obj.observation.z_obs = z_obs;
            obj.observation.grid.X = X;
            obj.observation.grid.Y = Y;
            obj.observation.grid.Z = Z;
        end
        
        
        % Method to define micro robots
        function obj = defineMicroRobots(obj, groupParams, defaultMagnetParams)
            % Combine robots from all groups into a flat structure
            allPositions = []; % Store all positions
            allMoments = [];   % Store all magnetic moments
            allMasses = [];    % Store all masses
            
            % Loop through each group to define micro robots
            for g = 1:numel(groupParams)
                % Extract parameters for the current group
                center = groupParams(g).center;
                gridSize = groupParams(g).gridSize;
                spacing = groupParams(g).spacing;
                
                % Create the grid of robot positions for this group
                [x_grid, y_grid] = meshgrid( ...
                    linspace(center(1) - (gridSize(2)-1)*spacing(1)/2, center(1) + (gridSize(2)-1)*spacing(1)/2, gridSize(2)), ...
                    linspace(center(2) - (gridSize(1)-1)*spacing(2)/2, center(2) + (gridSize(1)-1)*spacing(2)/2, gridSize(1)) ...
                );
                x_positions = reshape(x_grid, [], 1);
                y_positions = reshape(y_grid, [], 1);
                z_positions = zeros(size(x_positions)); % Assume z = 0 for all robots
                
                % Add positions to the combined list
                allPositions = [allPositions; [x_positions, y_positions, z_positions]];
                
                % Compute and add magnetic moments and masses
                V = prod(defaultMagnetParams.dimensions); % Volume of each robot
                m_moment_vector = defaultMagnetParams.M * V * defaultMagnetParams.orientation; % Magnetic moment
                allMoments = [allMoments; repmat(m_moment_vector, numel(x_positions), 1)];
                
                rho = defaultMagnetParams.rho; % Density of the material
                mass = rho * V; % Mass of each robot
                allMasses = [allMasses; repmat(mass, numel(x_positions), 1)];
            end
            
            % Store the combined flat list of micro robots
            obj.microRobots.positions = allPositions;
            obj.microRobots.moments = allMoments;
            obj.microRobots.masses = allMasses;
            obj.microRobots.defaultMagnetParams = defaultMagnetParams;
        end
        
        
        
        function obj = defineEquilibriumPoints(obj, times, positions)
            % Define equilibrium points at specific times and positions.
            % Inputs:
            %   times     - A vector of times at which equilibrium points are activated.
            %   positions - A Nx2 matrix of [x, y] positions for equilibrium points.

            if length(times) ~= size(positions, 1)
                error('Number of times must match the number of positions.');
            end

            % Initialize eq_points as an empty array of structs
            obj.eq_points = struct('time', {}, 'positions', {});

            % Group equilibrium points by unique time values
            uniqueTimes = unique(times); % Get unique time entries

            for tIdx = 1:length(uniqueTimes)
                currentTime = uniqueTimes(tIdx);
                % Find indices where the current time appears in the original time array
                matchIdx = find(times == currentTime);
                matchedPositions = positions(matchIdx, :);

                % If positions are identical, keep only one
                matchedPositions = unique(matchedPositions, 'rows');

                % Store in eq_points struct
                obj.eq_points(tIdx).time = currentTime;
                obj.eq_points(tIdx).positions = matchedPositions; % Store all positions active at this time
            end
        end

        
        
        
        
        % generate circular equilibrium points
        function [new_eq_positions, new_eq_times] = generateCircularEqPoints(obj, centerPoint, r1, N, s1, startTime)
            theta = linspace(0, 2*pi, N);  % Generate N angles
            circle_x = centerPoint(1) + r1 * cos(theta);
            circle_y = centerPoint(2) + r1 * sin(theta);
            
            new_eq_positions = [circle_x' circle_y']; % Convert to matrix
            new_eq_times = startTime + (1:N)' * s1;
        end
        
        
        % Method to generate points along a straight line
        function [new_eq_positions, new_eq_times] = generateLinearEqPoints(obj, startPoint, endPoint, N, s1, startTime)
            x_values = linspace(startPoint(1), endPoint(1), N);
            y_values = linspace(startPoint(2), endPoint(2), N);
            new_eq_positions = [x_values' y_values']; 
            new_eq_times = startTime + (1:N)' * s1;
        end
        
    end
end
