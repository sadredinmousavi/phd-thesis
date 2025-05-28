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

        
        
        
        % Method to save equilibrium data
        function saveEquilibriumData(obj, filename)
            timed_psai_array = obj.timed_psai_array;
            save(filename, 'timed_psai_array');
            fprintf('💾 Saved equilibrium data to %s\n', filename);
        end


        % Method to load equilibrium data
        function obj = loadEquilibriumData(obj, filename)
            if isfile(filename)
                loadedData = load(filename, 'timed_psai_array'); % Load stored variable
                obj.timed_psai_array = loadedData.timed_psai_array; % Assign to object property
                fprintf('🔄 Loaded equilibrium data from %s\n', filename);
            else
                fprintf('⚠ File not found: %s\n', filename);
            end
        end

        
        

        
        % Method to generate experimental inputs and save to file
        function generateExperimentalInputs(obj, filename)
            fileID = fopen(filename, 'w');
            for i = 1:length(obj.timed_psai_array)
                psai_values = obj.timed_psai_array(i).psai_array;
                rotation_degrees = rad2deg(psai_values);
                formatted_values = arrayfun(@(x) ...
                    sprintf('%s%.2f%s', x < 0 * '', x, x < 0 * ''), rotation_degrees, 'UniformOutput', false);

                formatted_command = sprintf('moveArray(port, [%s], 5, 180)\n', strjoin(formatted_values, ', '));
                fprintf(fileID, '%s', formatted_command);
            end
            fclose(fileID);
            disp("✅ Experimental inputs saved successfully.");
        end
        
        
        function plotStaticSnapshot(obj, psai_array, eq_positions)
            % plotStaticSnapshot Computes the force field using force_field_fast
            % for the provided psai_array and visualizes the result using streamslice.
            % Also displays psai angles, eq positions, and stability info below the plot.
            %
            % Inputs:
            %   psai_array - Magnet rotation angles. Can be manually provided or
            %                selected from obj.timed_psai_array.
            %   eq_positions - Equilibrium positions to analyze (Nx2 matrix)

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

            % Create the figure with subplots
            fig = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

            % Main plot area (top 70% of figure)
            ax = subplot(2, 1, 1);

            % Visualize force field using streamslice with cubic interpolation
            plot_field = streamslice(ax, X, Y, Fx, Fy, 'method', 'cubic');
            set(plot_field, 'Color', 'black', 'LineWidth', 1.2);
            hold(ax, 'on');

            % Overlay magnets
            mag_positions = obj.defs.epms.positions;
            plot(ax, mag_positions(:,1), mag_positions(:,2), 'r.', 'MarkerSize', 15);

            % Plot equilibrium points if provided
            if ~isempty(eq_positions)
                if size(eq_positions, 2) ~= 2
                    error('eq_positions must be Nx2 matrix');
                end
                plot(ax, eq_positions(:,1), eq_positions(:,2), 'rx', ...
                     'MarkerSize', 15, 'LineWidth', 2);
            end

            % Apply formatting and axis settings dynamically
            title(ax, '\textbf{Force Field}', 'interpreter', 'latex', 'fontsize', 14);
            xlabel(ax, 'x [m]', 'interpreter', 'latex');
            ylabel(ax, 'y [m]', 'interpreter', 'latex');

            % Calculate domain limits inline
            x_min = min(obj.observationSpace.grid.X(:));
            x_max = max(obj.observationSpace.grid.X(:));
            y_min = min(obj.observationSpace.grid.Y(:));
            y_max = max(obj.observationSpace.grid.Y(:));

            % Set axis limits based on observation space
            xlim(ax, [x_min, x_max]);
            ylim(ax, [y_min, y_max]);
            axis(ax, 'square');
            grid(ax, 'on');
            hold(ax, 'off');

            % Information display area (bottom 30% of figure)
            info_ax = subplot(2, 1, 2);
            axis(info_ax, 'off');

            % Create text content
            psai_text = sprintf('Magnet Angles (rad):\n');
            for i = 1:length(psai_array)
                psai_text = sprintf('%sMagnet %d: %.3f°\n', psai_text, i, rad2deg(psai_array(i)));
            end

            if ~isempty(eq_positions)
                eq_text = sprintf('\nEquilibrium Positions:\n');
                for i = 1:size(eq_positions,1)
                    eq_text = sprintf('%s(%.3f, %.3f)\n', eq_text, eq_positions(i,1), eq_positions(i,2));
                end

                stability_text = sprintf('\nStability Analysis:\n');
                for i = 1:size(eq_positions,1)
                    % Compute Hessian at equilibrium point
                    [~, H] = force_field_fast(eq_positions(i,1), eq_positions(i,2), psai_array);

                    % Compute eigenvalues and eigenvectors
                    [V, D] = eig(H);
                    eigenvals = diag(D);

                    stability_text = sprintf('%sPoint (%.2f, %.2f):\n', stability_text, ...
                                           eq_positions(i,1), eq_positions(i,2));
                    stability_text = sprintf('%s  Eigenvalues: %.3e, %.3e\n', stability_text, ...
                                   eigenvals(1), eigenvals(2));
                    stability_text = sprintf('%s  Eigenvectors:\n', stability_text);
                    stability_text = sprintf('%s    [%.3f; %.3f]\n', stability_text, V(1,1), V(2,1));
                    stability_text = sprintf('%s    [%.3f; %.3f]\n\n', stability_text, V(1,2), V(2,2));
                end
            else
                eq_text = '\nNo equilibrium points specified';
                stability_text = '\nNo equilibrium points to analyze';
            end

            % Combine all text
            full_text = [psai_text eq_text stability_text];

            % Display text
            text(0, 0.5, full_text, 'FontName', 'FixedWidth', 'FontSize', 10, ...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
                'Interpreter', 'none');

            % Adjust subplot positions for RTL layout
            set(info_ax, 'Position', [0.05 0.1 0.25 0.8]);  % Left panel (info)
            set(ax, 'Position', [0.35 0.1 0.6 0.8]);       % Right panel (plot)
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
            % createSimulationVideo Creates optimized video with force field visualization option
            % New Parameter:
            %   'ShowForceField' - Toggle force field display (default: false)

            % Parse inputs
            p = inputParser;
            addParameter(p, 'FrameRate', 10, @isnumeric);
            addParameter(p, 'ShowMagnets', true, @islogical);
            addParameter(p, 'FrameSkip', 1, @isnumeric);
            addParameter(p, 'ShowForceField', false, @islogical); % New parameter
            parse(p, varargin{:});

            % Verify data
            if isempty(obj.robot_trajectories)
                error('No trajectories available. Run simulation first.');
            end

            % Create figure with faster rendering
            f = figure('Position', [100 100 800 600], 'Color', 'w', 'Visible', 'off');
            ax = axes('Parent', f);
            hold(ax, 'on');

            % Set axis limits
            x_lim = [min(obj.defs.observation.x_space), max(obj.defs.observation.x_space)] + [-0.05 0.05];
            y_lim = [min(obj.defs.observation.y_space), max(obj.defs.observation.y_space)] + [-0.05 0.05];
            axis(ax, 'equal');
            xlim(ax, x_lim);
            ylim(ax, y_lim);

            % Initialize force field visualization (if enabled)
            if p.Results.ShowForceField
                [X, Y] = meshgrid(linspace(x_lim(1), x_lim(2), 20), linspace(y_lim(1), y_lim(2), 20));
                Fx = zeros(size(X));
                Fy = zeros(size(Y));
                force_plot = streamslice(ax, X, Y, Fx, Fy, 'method', 'cubic');
                set(force_plot, 'Color', 'black', 'LineWidth', 1.2);
            end

            % Plot magnets if enabled
            if p.Results.ShowMagnets
                mag_pos = obj.defs.epms.positions;
                plot(ax, mag_pos(:,1), mag_pos(:,2), 'r^', ...
                     'MarkerSize', 10, 'LineWidth', 2);
            end

            % Initialize robot plots
            robot_plots = gobjects(length(obj.robot_trajectories), 1);
            for i = 1:length(obj.robot_trajectories)
                robot_plots(i) = plot(ax, NaN, NaN, 'o', ...
                    'MarkerFaceColor', [0 0.447 0.741], ...
                    'MarkerEdgeColor', 'k', ...
                    'MarkerSize', 6);
            end

            % Video writer setup
            v = VideoWriter(filename, 'MPEG-4');
            v.FrameRate = p.Results.FrameRate;
            v.Quality = 90;
            open(v);

            % Main rendering loop
            frame_times = 0:(1/p.Results.FrameRate):obj.simulation_time(end);
            h = waitbar(0, 'Rendering video...');

            for k = 1:p.Results.FrameSkip:length(frame_times)
                current_time = frame_times(k);
                all_psai_times = [obj.timed_psai_array.time]; % Extract stored equilibrium times
                [~, psai_idx] = min(abs(all_psai_times - current_time)); % Find closest match
                % Retrieve the current psai configuration dynamically
                psai_array = obj.timed_psai_array(psai_idx).psai_array;
                
                [~, idx] = min(abs(obj.simulation_time - current_time));
                
                % Update robot positions
                for i = 1:length(obj.robot_trajectories)
                    set(robot_plots(i), ...
                        'XData', obj.robot_trajectories(i).x(idx), ...
                        'YData', obj.robot_trajectories(i).y(idx));
                end
                
                % Update force field dynamically
                if p.Results.ShowForceField
                    for xIdx = 1:size(X, 1)
                        for yIdx = 1:size(X, 2)
                            x_val = X(xIdx, yIdx);
                            y_val = Y(xIdx, yIdx);

                            % Compute force using the dynamically retrieved psai_array
                            F_total = force_field_fast(x_val, y_val, psai_array);
                            Fx(xIdx, yIdx) = F_total(1);
                            Fy(xIdx, yIdx) = F_total(2);
                        end
                    end

                    % Remove previous field & update visualization
                    delete(force_plot);  
                    force_plot = streamslice(ax, X, Y, Fx, Fy, 2);
                    set(force_plot, 'Color', 'black', 'LineWidth', 1.2);
                end

                % Update title
                title(ax, sprintf('Time = %.2f s', current_time));

                % Capture frame
                writeVideo(v, getframe(f));

                % Update progress
                waitbar(k/length(frame_times), h);
            end

            % Cleanup
            close(v);
            close(f);
            close(h);
            fprintf('Video saved: %s\n', filename);
        end

        
        
        
        

%         function createSimulationVideo(obj, filename, varargin)
%             % createSimulationVideo Creates real-time video with 0.1s resolution
%             % Inputs:
%             %   filename - Output video filename
%             % Optional:
%             %   'FrameRate' - Target frame rate (default: 10 for 0.1s steps)
%             %   'ShowMagnets' - Toggle magnet display (default: true)
% 
%             % Parse inputs with simplified options
%             p = inputParser;
%             addParameter(p, 'FrameRate', 10, @isnumeric); % 10 fps = 0.1s steps
%             addParameter(p, 'ShowMagnets', true, @islogical);
%             parse(p, varargin{:});
% 
%             % Verify data
%             if isempty(obj.robot_trajectories)
%                 error('No trajectories available. Run simulation first.');
%             end
% 
%             % Create video writer
%             v = VideoWriter(filename, 'MPEG-4');
%             v.FrameRate = p.Results.FrameRate;
%             v.Quality = 90;
%             open(v);
% 
%             % Create figure with optimized settings
%             f = figure('Position', [100 100 800 600], 'Color', 'w', 'Visible', 'off');
%             ax = axes('Parent', f, 'SortMethod', 'childorder');
%             hold(ax, 'on');
% 
%             % Set fixed axis limits
%             xlim(ax, [min(obj.defs.observation.x_space), max(obj.defs.observation.x_space)] + [-0.05 0.05]);
%             ylim(ax, [min(obj.defs.observation.y_space), max(obj.defs.observation.y_space)] + [-0.05 0.05]);
% 
%             % Plot magnets if enabled (red triangles)
%             if p.Results.ShowMagnets
%                 mag_pos = obj.defs.epms.positions;
%                 plot(ax, mag_pos(:,1), mag_pos(:,2), 'r^', 'MarkerSize', 10, 'LineWidth', 2);
%             end
% 
%             % Initialize robot plots (all blue, 0.005m diameter)
%             robot_color = [0 0.447 0.741]; % Standard MATLAB blue
%             robot_plots = gobjects(length(obj.robot_trajectories), 1);
% 
%             for i = 1:length(obj.robot_trajectories)
%                 robot_plots(i) = plot(ax, NaN, NaN, 'o', ...
%                     'MarkerFaceColor', robot_color, ...
%                     'MarkerEdgeColor', 'k', ...
%                     'MarkerSize', 6, ... % ~0.005m when figure is properly scaled
%                     'LineWidth', 0.5);
%             end
% 
%             % Configure axes
%             axis(ax, 'equal');
%             grid(ax, 'on');
%             xlabel(ax, 'X position [m]');
%             ylabel(ax, 'Y position [m]');
%             title(ax, 'Time = 0.0 s');
% 
%             % REAL-TIME RENDERING LOOP (0.1s resolution)
%             time_step = 0.3; % Fixed 0.1s resolution
%             current_time = 0;
%             k = 1; % Simulation index
%             last_frame_time = 0;
% 
%             while k <= length(obj.simulation_time)
%                 % Find next simulation point >= current_time
%                 while k <= length(obj.simulation_time) && obj.simulation_time(k) < current_time
%                     k = k + 1;
%                 end
% 
%                 if k > length(obj.simulation_time)
%                     break; % Reached end of simulation
%                 end
% 
%                 % Update visualization
%                 title(ax, sprintf('Time = %.1f s', current_time));
%                 for i = 1:length(obj.robot_trajectories)
%                     set(robot_plots(i), ...
%                         'XData', obj.robot_trajectories(i).x(k), ...
%                         'YData', obj.robot_trajectories(i).y(k));
%                 end
% 
%                 % Capture frame exactly at 0.1s intervals
%                 if current_time >= last_frame_time
%                     writeVideo(v, getframe(f));
%                     last_frame_time = last_frame_time + (1/p.Results.FrameRate);
%                 end
% 
%                 drawnow limitrate;
%                 current_time = current_time + time_step;
%             end
% 
%             % Close video file
%             close(v);
%             close(f);
%             fprintf('Video saved: %s (%.1f s duration, %.1f fps)\n', ...
%             filename, obj.simulation_time(end), p.Results.FrameRate);
%         end
        
       






























        % calcEquilibrium calculates an optimal psai_array for each equilibrium point.
        % It uses fmincon and a cost function based on force_field_fast.
        function obj = calcEquilibrium(obj, useSavedData, filename)
            if useSavedData && isfile(filename) % Check if file exists
                fprintf('✅ Loading precomputed equilibrium data...\n');
                obj = obj.loadEquilibriumData(filename);
                return;
            end
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

                % Step 1: First Optimization - Force Minimization
                % Define cost function that minimizes force at all target points simultaneously
                costFun = @(psai) obj.computeCost(psai, target_positions);

                % Improved optimization setup
                options = optimoptions('fmincon',...
                    'SpecifyObjectiveGradient', true,...
                    'Display', 'iter-detailed',...
                    'Algorithm', 'interior-point', ...  % Better for nonlinear constraints
                    'MaxIterations', 5000,...
                    'TolFun', 1e-9, ...    % Tighter tolerance
                    'TolCon', 1e-9, ...     % Constraint tolerance
                    'TolX', 1e-9,...
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
                [opt_psai, fval] = run(ms, problem, 5);  % 5 random starts

                % Perform optimization
                %[opt_psai, fval] = fmincon(costFun, psai0, [], [], [], [], lb, ub, [], options);
                
                fprintf('First Optimization Complete: Cost = %.4e\n', fval);
                
                
                % Step 2: Refinement - Eigenvalue Similarity Optimization
                costFun_zoomed = @(psai) obj.computeCost_zoomed(psai, target_positions);
                refinement_lb = max(opt_psai * 0.95, lb); % Lower bound: No more than -5% change
                refinement_ub = min(opt_psai * 1.05, ub); % Upper bound: No more than +5% change
                
                % Run fmincon with refined bounds
                options_zoomed = optimoptions('fmincon', ...
                    'Display', 'iter-detailed', ...
                    'Algorithm', 'interior-point', ...
                    'MaxIterations', 5000, ...
                    'TolFun', 1e-9, ...
                    'TolX', 1e-9);

                opt_psai_zoomed = fmincon(costFun_zoomed, opt_psai, [], [], [], [], refinement_lb, refinement_ub, [], options_zoomed);
        
        
                % Store the optimized psai_array for this equilibrium time
                obj.timed_psai_array(idx).time = target_time;
                obj.timed_psai_array(idx).psai_array = opt_psai_zoomed; % Store one psai solution for all points

                fprintf('Equilibrium at time %.2f, positions:\n', target_time);
                disp(target_positions);
                fprintf('Optimized psai_array:\n');
                disp(opt_psai);
                fprintf('Cost: %.4e\n', fval);
            end

            fprintf('Equilibrium calculation complete.\n');
            % Save results for future use
            obj.saveEquilibriumData(filename);
            fprintf('💾 Equilibrium data saved!\n');
        end
        
      
        
        
        
        
        
        
        
        
        function [cost, grad] = computeCost(obj, psai, target_positions)
            cost = 0;
            grad = zeros(size(psai)); % grad should be 1x8
            for i = 1:size(target_positions,1)
                [F, Grad_psai] = force_gradients_fast(target_positions(i,1), target_positions(i,2), psai);
                F_xy = F(1:2); % Make it row vector (1x2)
                J = Grad_psai(:,1:2); % Jacobian (8x2)

                % Calculate the norm of the force vector
                norm_F = sqrt(sum(F_xy.^2));

                % Cost calculation: sum of force magnitudes
                cost = cost + norm_F;

                % Gradient calculation using the chain rule:
                % derivative of sqrt(s) with respect to s is 1/(2*sqrt(s)) and ds/dF_xy = 2*F_xy,
                % so the derivative simplifies to F_xy/norm_F.
                if norm_F > 1e-12  % Avoid division by zero for zero force
                    grad = grad + (F_xy / norm_F) * J';
                end
            end
        end
        
        function cost = computeCost_zoomed(obj, psai, target_positions)
            cost = 0;
            weight_eig = 10; % Scaling factor for eigenvalue similarity

            for i = 1:size(target_positions,1)
                % Compute Hessian for the current position
                [~, H] = force_field_fast(target_positions(i,1), target_positions(i,2), psai);
                eigen_vals = eig(H);  % Extract eigenvalues

                if min(eigen_vals) > 1e-12
                    ratio = max(eigen_vals) / min(eigen_vals); % Ratio between largest and smallest eigenvalue
                    cost = cost + weight_eig * (ratio - 1)^2; % Penalize deviation from 1 (perfectly equal)
                else
                    % If eigenvalues are too small, just penalize based on squared difference
                    cost = cost + weight_eig * (eigen_vals(1) - eigen_vals(2))^2;
                end
            end
        end

        
%         function [cost, grad] = computeCost(obj, psai, target_positions)
%             cost = 0;
%             grad = zeros(size(psai));  % grad should be 1x8
%             weight_eig = 10000;
% 
%             for i = 1:size(target_positions,1)
%                 %%%--- Existing Force and Gradient Calculation ---%%%
%                 [F, Grad_psai] = force_gradients_fast(target_positions(i,1), target_positions(i,2), psai);
%                 F_xy = F(1:2);  % 1x2 vector (force in x and y)
%                 % Assume Grad_psai(:,1:2) contains the Jacobian relating F_xy to psai
%                 J = Grad_psai(:,1:2);
% 
%                 norm_F = sqrt(sum(F_xy.^2));
%                 cost = cost + norm_F;
%                 if norm_F > 1e-12
%                     % Chain rule calculation (derivative of sqrt(F^2) gives F/||F||)
%                     grad = grad + (F_xy / norm_F) * J';
%                 end
% 
%                 %%%--- New Hessian Eigenvalue Penalty ---%%%
%                 % Use the force_field_fast function (which returns Hessian)
%                 [~, H] = force_field_fast(target_positions(i,1), target_positions(i,2), psai);
%                 eigen_vals = eig(H);  % For the 2×2 Hessian, this returns a 2×1 vector.
%                 % Compute the squared difference of the eigenvalues.
%                 penalty = weight_eig * (eigen_vals(1) - eigen_vals(2))^2;
%                 cost = cost + penalty;
% 
%                 % Note: Deriving the gradient contribution of this penalty term is non‐trivial.
%                 % You have two options:
%                 %    (a) If weight_eig is chosen very small, you might assume its gradient is negligible.
%                 %    (b) Alternatively, derive the sensitivity of the eigenvalues (using eigenvalue
%                 %        perturbation theory) and propagate through psai.
%                 % For now, we are not adding an analytical gradient for the penalty term.
%             end
% 
% 
% %             cost = 0;
% %             grad = zeros(size(psai));
% %             force_scale = 1e-6;    % Typical force magnitude (adjust empirically)
% %             min_eig = 1e-12;        % Minimum eigenvalue to avoid division by zero
% %             epsilon = 1e-7;        % Smoothing term for weight transition
% % 
% %             for i = 1:size(target_positions, 1)
% %                 % --- Force Term ---
% %                 [F, Grad_psai] = force_gradients_fast(target_positions(i,1), target_positions(i,2), psai);
% %                 F_xy = F(1:2);
% %                 J = Grad_psai(:,1:2);  % Jacobian (8×2)
% % 
% %                 % Normalized force
% %                 norm_F = norm(F_xy) / force_scale;
% %                 cost_force = norm_F;
% % 
% %                 % Force gradient (chain rule)
% %                 if norm_F > 1e-12
% %                     grad_force = (F_xy / (norm_F * force_scale)) * J';
% %                 else
% %                     grad_force = 0;
% %                 end
% % 
% %                 % --- Eigenvalue Ratio Term ---
% %                 [~, H] = force_field_fast(target_positions(i,1), target_positions(i,2), psai);
% %                 eigen_vals = eig(H);
% %                 if numel(eigen_vals) >= 2 && all(abs(eigen_vals) > min_eig)
% %                     % Symmetric relative difference: (λ1 - λ2) / (0.5*(λ1 + λ2))
% %                     eig_avg = 0.5 * (eigen_vals(1) + eigen_vals(2));
% %                     rel_diff = (eigen_vals(1) - eigen_vals(2)) / (eig_avg + min_eig);
% %                     cost_eig = rel_diff^2;
% % 
% %                     % Gradient of eigenvalue term (simplified)
% % %                     [V, ~] = eig(H);
% % %                     dH_dpsai = Grad_psai(:,3:end);  % Columns 3-4 for Hessian gradients
% % %                     drel_diff_dH = (V(:,1)*V(:,1)' - V(:,2)*V(:,2)') / (eig_avg + min_eig);
% % %                     grad_eig = 2 * rel_diff * drel_diff_dH(:)' * dH_dpsai;
% %                 else
% %                     cost_eig = 0;
% % %                     grad_eig = 0;
% %                 end
% % 
% %                 % --- Adaptive Weighting ---
% %                 % Weight increases as force decreases: weight_eig = 1 / (norm_F + ε)
% %                 weight_eig = 1 / (norm_F + epsilon);
% %                 cost = cost + cost_force + weight_eig * cost_eig;
% % %                 grad = grad + grad_force + weight_eig * grad_eig;
% %                 grad = grad + grad_force;
% %             end
%         end
        
        
        
    end
end
