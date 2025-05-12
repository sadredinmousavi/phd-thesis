classdef Simulation
    properties
        defs                % Instance of the Definitions class
        timed_psai_array    % Structure array with fields: time, psai_array (one per eq_point)
        observationSpace    % Read from defs.observation
        
        ode_results       % Stores all ODE simulation results
        robot_trajectories % Processed trajectories for each robot
        simulation_time   % Time vector used in simulation
    end
    methods
        % Constructor: stores definitions and observation space.
        function obj = Simulation(defs)
            obj.defs = defs;
            obj.timed_psai_array = struct('time', {}, 'psai_array', {});
            obj.observationSpace = defs.observation;
        end

        % calcEquilibrium calculates an optimal psai_array for each equilibrium point.
        % It uses fmincon and a cost function based on force_field_fast.
        function obj = calcEquilibrium(obj)
            % Retrieve equilibrium times and positions from Definitions
            eq_times = [obj.defs.eq_points.time];
            eq_positions = {obj.defs.eq_points.positions}; % Cell array for multiple positions per time
            numMagnets = obj.defs.epms.inputs.numMagnets;

            fprintf('Calculating psai_array for each equilibrium configuration...\n');

            for idx = 1:length(eq_times)
                target_time = eq_times(idx);      % Current equilibrium time
                target_positions = eq_positions{idx}; % Multiple [x, y] positions at this time

                % Optimization bounds for psai values
                lb = -90 * pi/180 * ones(1, numMagnets);
                ub = +70 * pi/180 * ones(1, numMagnets);
                psai0 = 0.2 * lb; % Small initial guess

                % Define cost function that minimizes force at all target points simultaneously
                costFun = @(psai) obj.computeCost(psai, target_positions);


                

                % Improved optimization setup
                options = optimoptions('fmincon',...
                    'SpecifyObjectiveGradient', true,...
                    'Display', 'final-detailed',...
                    'Algorithm', 'interior-point', ...  % Better for nonlinear constraints
                    'MaxIterations', 5000,...
                    'TolFun', 1e-16, ...    % Tighter tolerance
                    'TolCon', 1e-16, ...     % Constraint tolerance
                    'TolX', 1e-16,...
                    'FiniteDifferenceType', 'central', ... % More accurate derivatives
                    'FiniteDifferenceStepSize', 1e-3);   % Better gradient estimation

                % Multi-start optimization (helps avoid local minima)
                ms = MultiStart('UseParallel', true);
                problem = createOptimProblem('fmincon',...
                    'objective', costFun,...
                    'x0', psai0,...
                    'lb', lb,...
                    'ub', ub,...
                    'options', options);
                [opt_psai, fval] = run(ms, problem, 20);  % 5 random starts

                % Perform optimization
                %[opt_psai, fval] = fmincon(costFun, psai0, [], [], [], [], lb, ub, [], options);

                % Store the optimized psai_array for this equilibrium time
                obj.timed_psai_array(idx).time = target_time;
                obj.timed_psai_array(idx).psai_array = opt_psai; % Store one psai solution for all points

                fprintf('Equilibrium at time %.2f, positions:\n', target_time);
                disp(target_positions);
                fprintf('Optimized psai_array:\n');
                disp(opt_psai);
                fprintf('Cost: %.4e\n', fval);
            end

            fprintf('Equilibrium calculation complete.\n');
        end

        
        
        function plotStaticSnapshot(obj, psai_array, eq_positions)
            % plotStaticSnapshot Computes the force field using force_field_fast
            % for the provided psai_array and visualizes the result using streamslice.
            %
            % Inputs:
            %   psai_array - Magnet rotation angles. Can be manually provided or
            %                selected from obj.timed_psai_array.

            if nargin < 3
                eq_positions = []; % Make it optional
            end
            
            % Retrieve observation grid from Definitions
            X = obj.observationSpace.grid.X;
            Y = obj.observationSpace.grid.Y;
            Z = obj.observationSpace.grid.Z; % Should be constant (e.g., Z = 0)

            % Initialize force field arrays
            Fx = zeros(size(X));
            Fy = zeros(size(Y));
            Fz = zeros(size(Z));  % Optional, should ideally be zero

            % Compute force field over the workspace
            for xIdx = 1:size(X, 1)
                for yIdx = 1:size(X, 2)
                    x_val = X(xIdx, yIdx);
                    y_val = Y(xIdx, yIdx);

                    % Compute force using the hardcoded function
                    F_total = force_field_fast(x_val, y_val, psai_array);

                    % Store force components
                    Fx(xIdx, yIdx) = F_total(1);
                    Fy(xIdx, yIdx) = F_total(2);
                    Fz(xIdx, yIdx) = F_total(3); % Likely zero
                end
            end

            % Create the figure and set its position manually
            fig = figure;
            screenSize = get(0, 'ScreenSize'); 
            fig.Position = [screenSize(1), screenSize(2), screenSize(3)*0.9, screenSize(4)*0.9]; 

            % Visualize force field using streamslice with cubic interpolation
            plot_field = streamslice(X, Y, Fx, Fy, 'method', 'cubic');
            set(plot_field, 'Color', 'black', 'LineWidth', 1.2);
            hold on;

            % Overlay magnets
            mag_positions = obj.defs.epms.positions;
            plot(mag_positions(:,1), mag_positions(:,2), 'r.', 'MarkerSize', 15);

            % Plot the equilibrium point only if it is provided
            % Plot equilibrium points if provided
            if ~isempty(eq_positions)
                if size(eq_positions, 2) ~= 2
                    error('eq_positions must be Nx2 matrix');
                end
                plot(eq_positions(:,1), eq_positions(:,2), 'rx', ...
                     'MarkerSize', 15, 'LineWidth', 2);
            end

            % Apply formatting and axis settings dynamically
            title('\textbf{Force Field}', 'interpreter', 'latex', 'fontsize', 14);
            xlabel('x [m]', 'interpreter', 'latex');
            ylabel('y [m]', 'interpreter', 'latex');

            % Calculate domain limits inline
            x_min = min(obj.observationSpace.grid.X(:));
            x_max = max(obj.observationSpace.grid.X(:));
            y_min = min(obj.observationSpace.grid.Y(:));
            y_max = max(obj.observationSpace.grid.Y(:));

            % Set axis limits based on observation space
            xlim([x_min, x_max]);
            ylim([y_min, y_max]);
            axis square;

            % Finalize plot
            hold off;
            grid on;
        end

        
        
        
        
        
        
        function obj = runODESimulation(obj, tspan, initial_conditions, options)
            % runODESimulation Solves ODE for microrobot motion under magnetic forces
            % Inputs:
            %   tspan - [t0 tf] or vector of time points
            %   initial_conditions - Nx2 matrix of [x,y] positions for each robot
            %   options - odeset options (optional)

            if nargin < 4
                options = odeset('RelTol',1e-6,'AbsTol',1e-9);
            end

            % Get number of robots
            num_robots = size(obj.defs.microRobots.positions, 1);
            
            % Verify initial conditions match number of robots
            if size(initial_conditions,1) ~= num_robots
                error('Initial conditions must have %d rows (one per robot)', num_robots);
            end
            
            % Verify initial conditions now expect [x,y,vx,vy] for each robot
            if size(initial_conditions, 2) ~= 4
                error('Initial conditions must be Nx4 matrix: [x,y,vx,vy] for each robot');
            end

            % Flatten initial conditions [x1,y1,vx1,vy1, x2,y2,vx2,vy2, ...]
            y0 = reshape(initial_conditions', [], 1);

            % Get default magnet angles (or use current psai_array if available)
            if ~isempty(obj.timed_psai_array)
                psai_array = obj.timed_psai_array(1).psai_array; % Use first configuration
            else
                psai_array = zeros(1, obj.defs.epms.inputs.numMagnets);
            end
    
            % Get equilibrium points data
            eq_times = [obj.defs.eq_points.time];
            [sorted_times, time_idx] = sort(eq_times);
            psai_configs = {obj.timed_psai_array(time_idx).psai_array};

            % Define ODE function with time-dependent psai
            function dydt = ode_func(t, y)
                % Find current psai configuration
                current_psai = zeros(1, obj.defs.epms.inputs.numMagnets); % Default
                for k = length(sorted_times):-1:1
                    if t >= sorted_times(k)
                        current_psai = psai_configs{k};
                        break;
                    end
                end

                % Calculate derivatives
                dydt = obj.systemDynamics(t, y, current_psai);
            end

            % Solve ODE
            [t, y] = ode45(@ode_func, tspan, y0, options);

            % Store results
            obj.ode_results.t = t;
            obj.ode_results.y = y;
            obj.simulation_time = t;

            % Process trajectories into more usable format
            obj = obj.processTrajectories();
        end

        
        
        
        function [cost, grad] = computeCost(obj, psai, target_positions)
            cost = 0;
            grad = zeros(size(psai)); % grad should be 1x8
            for i = 1:size(target_positions,1)
                [F, Grad_psai] = force_gradients_fast(target_positions(i,1), target_positions(i,2), psai);
                F_xy = F(1:2); % Make it row vector (1x2)
                J = Grad_psai(:,1:2); % Jacobian (8x2)

                % Cost calculation (sum of squared forces)
                cost = cost + sum(F_xy.^2);

                % Gradient calculation (1x8) = (1x2) * (2x8)
                grad = grad + 2 * F_xy * J';
            end
        end
        
        
        
        function dydt = systemDynamics(obj, t, y, psai_array)
    % systemDynamics Complete dynamics including:
    % - Magnetic forces (external field)
    % - Dipole-dipole interactions
    % - Drag forces
    
    % Get parameters from definitions
    mu_0 = obj.defs.epms.inputs.mu_0;
    m = obj.defs.microRobots.defaultMagnetParams.M * ...
        prod(obj.defs.microRobots.defaultMagnetParams.dimensions);
    diameter = obj.defs.microRobots.defaultMagnetParams.dimensions(1);
    Mu = obj.defs.microRobots.defaultMagnetParams.Mu;
    rho = obj.defs.microRobots.defaultMagnetParams.rho;
    volume = prod(obj.defs.microRobots.defaultMagnetParams.dimensions);
    mass = rho * volume;
    
    % Reshape state vector: [x1;y1;vx1;vy1; x2;y2;vx2;vy2; ...]
    state_matrix = reshape(y, 4, [])';
    num_robots = size(state_matrix, 1);
    
    % Initialize derivatives
    dydt = zeros(size(y));
    
    % Calculate all positions and velocities
    positions = state_matrix(:,1:2);
    velocities = state_matrix(:,3:4);
    
    % Calculate interaction forces
    F_inter = zeros(num_robots, 2);
    for i = 1:num_robots
        for j = (i+1):num_robots
            r_vec = positions(i,:) - positions(j,:);
            r = max(norm(r_vec), 1e-6); % Avoid division by zero
            r_hat = r_vec/r;
            
            % Dipole-dipole force (simplified 2D)
            force_mag = (3*mu_0*m^2)/(4*pi*r^4);
            F_ij = force_mag * r_hat;
            
            F_inter(i,:) = F_inter(i,:) + F_ij;
            F_inter(j,:) = F_inter(j,:) - F_ij;
        end
    end
    
    % Calculate derivatives for each robot
    for i = 1:num_robots
        idx = (i-1)*4 + (1:4);
        pos = positions(i,:);
        vel = velocities(i,:);
        
        % 1. External magnetic force
        F_magnetic = force_field_fast(pos(1), pos(2), psai_array);
        
        % 2. Drag force (Stokes)
        drag_coeff = Mu * (16/3) * diameter;
        F_drag = -drag_coeff * vel';
        
        % 3. Sum all forces
        F_total = F_magnetic(1:2)' + F_drag + F_inter(i,:)';
        
        % State derivatives
        dydt(idx(1)) = vel(1);          % dx/dt = vx
        dydt(idx(2)) = vel(2);          % dy/dt = vy
        dydt(idx(3)) = F_total(1)/mass; % dvx/dt = Fx/m
        dydt(idx(4)) = F_total(2)/mass; % dvy/dt = Fy/m
    end
        end
        
        
        

        
        
        
        
        
        function obj = processTrajectories(obj)
    % PROCESSEDTRAJECTORIES Converts ODE results to position trajectories
    % Handles state vector format: [x1,y1,vx1,vy1, x2,y2,vx2,vy2, ...]
    
    % Check for results
    if isempty(obj.ode_results)
        error('No ODE results available. Run runODESimulation first.');
    end
    
    % Get basic parameters
    y = obj.ode_results.y;
    num_robots = size(obj.defs.microRobots.positions, 1);
    
    % Initialize trajectory storage
    obj.robot_trajectories = struct('x', [], 'y', []);
    
    % Extract positions (skip velocities)
    for i = 1:num_robots
        % Position indices in state vector
        x_idx = (i-1)*4 + 1;
        y_idx = (i-1)*4 + 2;
        
        % Store trajectories
        obj.robot_trajectories(i).x = y(:, x_idx);
        obj.robot_trajectories(i).y = y(:, y_idx);
        
        % Debug check (optional)
        if any(isnan(obj.robot_trajectories(i).x)) || any(isnan(obj.robot_trajectories(i).y))
            warning('NaN detected in robot %d trajectory', i);
        end
    end
    
    % Verify time alignment
    if length(obj.simulation_time) ~= size(y, 1)
        error('Time vector length mismatch with state data');
    end
        end
        
        
        
        

     function createSimulationVideo(obj, filename, varargin)
    % createSimulationVideo Creates real-time video with 0.1s resolution
    % Inputs:
    %   filename - Output video filename
    % Optional:
    %   'FrameRate' - Target frame rate (default: 10 for 0.1s steps)
    %   'ShowMagnets' - Toggle magnet display (default: true)

    % Parse inputs with simplified options
    p = inputParser;
    addParameter(p, 'FrameRate', 10, @isnumeric); % 10 fps = 0.1s steps
    addParameter(p, 'ShowMagnets', true, @islogical);
    parse(p, varargin{:});

    % Verify data
    if isempty(obj.robot_trajectories)
        error('No trajectories available. Run simulation first.');
    end

    % Create video writer
    v = VideoWriter(filename, 'MPEG-4');
    v.FrameRate = p.Results.FrameRate;
    v.Quality = 90;
    open(v);

    % Create figure with optimized settings
    f = figure('Position', [100 100 800 600], 'Color', 'w', 'Visible', 'off');
    ax = axes('Parent', f, 'SortMethod', 'childorder');
    hold(ax, 'on');
    
    % Set fixed axis limits
    xlim(ax, [min(obj.defs.observation.x_space), max(obj.defs.observation.x_space)] + [-0.05 0.05]);
    ylim(ax, [min(obj.defs.observation.y_space), max(obj.defs.observation.y_space)] + [-0.05 0.05]);

    % Plot magnets if enabled (red triangles)
    if p.Results.ShowMagnets
        mag_pos = obj.defs.epms.positions;
        plot(ax, mag_pos(:,1), mag_pos(:,2), 'r^', 'MarkerSize', 10, 'LineWidth', 2);
    end

    % Initialize robot plots (all blue, 0.005m diameter)
    robot_color = [0 0.447 0.741]; % Standard MATLAB blue
    robot_plots = gobjects(length(obj.robot_trajectories), 1);
    
    for i = 1:length(obj.robot_trajectories)
        robot_plots(i) = plot(ax, NaN, NaN, 'o', ...
            'MarkerFaceColor', robot_color, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 6, ... % ~0.005m when figure is properly scaled
            'LineWidth', 0.5);
    end

    % Configure axes
    axis(ax, 'equal');
    grid(ax, 'on');
    xlabel(ax, 'X position [m]');
    ylabel(ax, 'Y position [m]');
    title(ax, 'Time = 0.0 s');

    % REAL-TIME RENDERING LOOP (0.1s resolution)
    time_step = 0.1; % Fixed 0.1s resolution
    current_time = 0;
    k = 1; % Simulation index
    last_frame_time = 0;
    
    while k <= length(obj.simulation_time)
        % Find next simulation point >= current_time
        while k <= length(obj.simulation_time) && obj.simulation_time(k) < current_time
            k = k + 1;
        end
        
        if k > length(obj.simulation_time)
            break; % Reached end of simulation
        end

        % Update visualization
        title(ax, sprintf('Time = %.1f s', current_time));
        for i = 1:length(obj.robot_trajectories)
            set(robot_plots(i), ...
                'XData', obj.robot_trajectories(i).x(k), ...
                'YData', obj.robot_trajectories(i).y(k));
        end

        % Capture frame exactly at 0.1s intervals
        if current_time >= last_frame_time
            writeVideo(v, getframe(f));
            last_frame_time = last_frame_time + (1/p.Results.FrameRate);
        end
        
        drawnow limitrate;
        current_time = current_time + time_step;
    end

    % Close video file
    close(v);
    close(f);
    fprintf('Video saved: %s (%.1f s duration, %.1f fps)\n', ...
        filename, obj.simulation_time(end), p.Results.FrameRate);
end
        
       
        
        
        
        % runSimulation runs the simulation over the total duration. At the correct time,
        % it “switches” to the corresponding psai_array.
        function runSimulation(obj, simulation_time)
            % Retrieve equilibrium times and sort if necessary.
            eq_times = [obj.timed_psai_array.time];
            [eq_times, idx] = sort(eq_times);
            eq_positions = obj.defs.eq_points.position(idx, :);
            eq_psai_arrays = {obj.timed_psai_array(idx).psai_array};

            current_time = 0;
            dt = 0.01;
            eq_idx = 1;
            
            fprintf('Starting simulation...\n');
            while current_time < simulation_time
                % Switch to the next equilibrium (i.e. psai_array) if the time is reached.
                if eq_idx <= length(eq_times) && current_time >= eq_times(eq_idx)
                    fprintf('At time %.2f, switching to equilibrium point (%.2f, %.2f) with psai_array:\n', ...
                        current_time, eq_positions(eq_idx, 1), eq_positions(eq_idx, 2));
                    disp(eq_psai_arrays{eq_idx});
                    % (Here, you would update your magnet model with eq_psai_arrays{eq_idx})
                    eq_idx = eq_idx + 1;
                end
                
                % Advance simulation time.
                current_time = current_time + dt;
            end
            fprintf('Simulation completed.\n');
        end
    end
end
