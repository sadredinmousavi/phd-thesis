classdef PermanentMagnetField
    %PERMANENTMAGNETFIELD Calculates the dipole magnetic field from an array of
    %permanent magnets using the dipole-approximation.
    %   Magnet data is taken from a definitions object (defs). It assumes defs.epms contains:
    %       positions    - an N×3 array of magnet positions.
    %       orientations - an N×3 array of effective dipole moment vectors.
    %       psai         - a 1×N array of rotation angles (radians).
    %       inputs.mu0   - permeability of free space.
    %       inputs.Br    - remanence (Tesla).

    properties
        mu_0      % Permeability (read from defs)
        magnets  % Array of magnet structures with fields: position, orientation, Br, psi.
        defs
    end
    
    methods
        function obj = PermanentMagnetField(defs)
            %CONSTRUCTOR Create a PermanentMagnetField object using definitions.
            %   defs should be an instance of your Definitions class whose epms field
            %   now contains:
            %     positions    - positions of the magnets.
            %     orientations - effective dipole moment vectors.
            %     psai         - rotation angles (psi) for each magnet.
            %     inputs.mu0   - permeability.
            %     inputs.Br    - remanence (Tesla).
            if nargin > 0
                obj.defs = defs;
                % Read mu0 from definitions
                obj.mu_0 = defs.epms.inputs.mu_0;
                
                % Get magnet data from definitions
                positions    = defs.epms.positions;          % N×3 array (positions)
                orientations = defs.epms.orientations;       % N×3 array (effective moment vectors)
                psi          = defs.epms.psai;               % 1×N vector (rotation angles)
                Br           = defs.epms.inputs.Br;          % remanence (assumed common to all magnets)
                
                % Build the magnets structure array
                N = size(positions,1);
                magnetsStruct = repmat(struct('position', [], 'orientation', [], 'Br', [], 'psi', []), N, 1);
                for i = 1:N
                    magnetsStruct(i).position = positions(i,:);
                    magnetsStruct(i).orientation = orientations(i,:);
                    magnetsStruct(i).Br = Br;  % if Br varies per magnet, adjust accordingly
                    magnetsStruct(i).psi = psi(i);
                    magnetsStruct(i).V = defs.epms.inputs.V;
                end
                obj.magnets = magnetsStruct;
            else
                obj.magnets = [];
            end
        end
        
        
        
        function F_total = calculateForceOnMicrorobotInPlane(obj, P_robot, m_robot)
            % CALCULATEFORCEONMICROROBOTINPLANE computes the net in-plane (XY)
            % force on a microrobot placed at position P_robot with dipole moment
            % m_robot. This force is computed by summing contributions from all
            % the source magnets (obj.magnets) using a symbolic-derived formula:
            %
            %    f_m = (3*mu0/(4*pi))*(m_source*m_robot/(r^4))*cos(psai)*r_hat_xy
            %
            % where:
            %   - m_source and m_robot are the norms of the source and robot dipole moments,
            %   - r is the distance between the source and the robot,
            %   - psai is a rotation parameter from the magnet,
            %   - r_hat_xy is the unit vector in the XY plane.
            %
            % The output F_total is a 1×3 force vector [Fx, Fy, 0].

            F_total = [0, 0, 0];  % Initialize net force as zero.

            for i = 1:numel(obj.magnets)
                mag = obj.magnets(i);
                P_source = mag.position; % 1×3 position of this source magnet

                % Compute the effective source magnet moment vector.
                % Include dimensions if available.
                if isfield(mag, 'dimensions') && ~isempty(mag.dimensions)
                    m_source_vec = mag.Br / obj.mu_0 * prod(mag.dimensions) * mag.orientation;
                else
                    m_source_vec = mag.Br / obj.mu_0 * mag.orientation;
                end

                % Compute scalar norms of the moments.
                m_norm_source = norm(m_source_vec);
                m_norm_robot  = norm(m_robot);  % Ensure m_robot is 1×3.

                % Compute relative vector and its norm.
                r = P_robot - P_source;  % 1×3 vector from magnet to robot.
                r_norm = norm(r);
                if r_norm == 0
                    continue;  % Skip singularity if robot coincides with a magnet.
                end

                % Compute unit vector and then its XY projection.
                r_hat = r / r_norm;
                r_hat_xy = r_hat(1:2);  % The force will be computed only in XY.

                % Get rotation parameter psai for this magnet.
                if isfield(mag, 'psai')
                    psai_i = mag.psai;
                else
                    psai_i = 0;
                end

                % Compute the scalar factor according to the modified formula.
                scalar_term = (3 * obj.mu_0 / (4 * pi)) * (m_norm_source * m_norm_robot / (r_norm^4)) * cos(psai_i);
                f_i_xy = scalar_term * r_hat_xy;  % This is a 1x2 force vector in XY.
                f_i = [f_i_xy, 0];               % Extend to full 3D ([Fx, Fy, 0]).

                % Sum the contributions from each magnet.
                F_total = F_total + f_i;
            end
        end
        
        
        
        function H = calculateHessianForMicrorobotInPlane(obj, P_robot, m_robot)
            % CALCULATEHESSIANFORMICROROBOTINPLANE computes the Hessian matrix of the
            % in-plane (XY) force field at the microrobot's position.
            %
            % The Hessian matrix is defined as:
            %   H = [dFx/dx  dFx/dy;
            %        dFy/dx  dFy/dy]
            %
            % Inputs:
            %   P_robot : 1×3 position vector of the microrobot [x,y,z]
            %   m_robot : 1×3 dipole moment vector of the microrobot [mx,my,mz]
            %
            % Output:
            %   H : 2×2 Hessian matrix (N/m)

            % Initialize Hessian components
            H = zeros(2, 2);  % [Hxx, Hxy; Hyx, Hyy]

            for i = 1:numel(obj.magnets)
                mag = obj.magnets(i);
                P_source = mag.position;  % 1×3 source magnet position

                % Compute source dipole moment
                if isfield(mag, 'dimensions') && ~isempty(mag.dimensions)
                    m_source_vec = mag.Br / obj.mu_0 * prod(mag.dimensions) * mag.orientation;
                else
                    m_source_vec = mag.Br / obj.mu_0 * mag.orientation;
                end

                % Scalar dipole moments
                m_norm_source = norm(m_source_vec);
                m_norm_robot = norm(m_robot);

                % Relative position vector (XY only)
                r = P_robot(1:2) - P_source(1:2);  % [dx, dy]
                r_norm = norm(r);

                if r_norm == 0
                    continue;  % Avoid division by zero
                end

                % Precompute common terms
                dx = r(1);
                dy = r(2);
                r_sq = r_norm^2;
                r_5 = r_norm^5;
                r_7 = r_norm^7;

                % Get rotation parameter
                psai_i = mag.psai;  % Assume psai field exists

                % Compute coefficient (3μ₀/4π) * m_source * m_robot * cos(psai)
                C = (3 * obj.mu_0 / (4 * pi)) * (m_norm_source * m_norm_robot) * cos(psai_i);

                % Hessian contributions for this magnet
                Hxx_i = C * (1/r_5 - 5*dx^2/r_7);
                Hxy_i = C * (-5*dx*dy/r_7);
                Hyy_i = C * (1/r_5 - 5*dy^2/r_7);

                % Accumulate contributions
                H(1,1) = H(1,1) + Hxx_i;  % dFx/dx
                H(1,2) = H(1,2) + Hxy_i;  % dFx/dy
                H(2,1) = H(2,1) + Hxy_i;  % dFy/dx (symmetry)
                H(2,2) = H(2,2) + Hyy_i;  % dFy/dy
            end
        end
        
        
        
        
        
        function makeForceFunctionHradCode(obj)
            % MAKEFORCEFUNCTIONHARDCODE Generates optimized force+Hessian computation function
            % 
            % This method creates a MATLAB function 'force_field_fast.m' that computes:
            % - Net in-plane magnetic force [Fx, Fy, 0] on a microrobot
            % - Hessian matrix of the force field at the robot position
            %
            % The function unrolls all magnet calculations for maximum speed
            % --------------------------------------------------------------

            % ======================= SETUP ================================
            % Read permanent magnet and robot parameters from object definitions

            % Permanent magnet array parameters (epms)
            % Expected fields:
            % - positions: Nx3 array of magnet positions [x,y,z] (meters)
            % - psai: Nx1 array of rotation angles (radians)
            % - Br: Magnet remanence (A/m) - scalar (all magnets same)
            % - dimensions: 1x3 magnet dimensions [a,b,c] (meters)
            epms = obj.defs.epms;
            %mag = obj.magnets(i);

            % Microrobot parameters (pm_robot)
            % Expected fields:
            % - M: Robot magnetization (A/m)
            % - dimensions: 1x3 robot magnet dimensions [a,b,c] (meters)
            pm_robot = obj.defs.microRobots.defaultMagnetParams;

            % ================== PHYSICAL CONSTANTS ========================
            % Compute magnetic moments using magnet volumes and material properties

            % Vacuum permeability (T·m/A)
            mu0 = obj.mu_0;  % Read from parent object

            % Robot dipole moment (A·m²)
            V_robot = prod(pm_robot.dimensions);    % Volume (m³)
            m_robot = pm_robot.M * V_robot;         % Magnetic moment (A·m²)

            % Source magnets dipole moment (same for all magnets)
            V_source = epms.inputs.V;                      % Magnet volume (m³)
            m_source = (epms.inputs.Br / mu0) * V_source;  % Magnetic moment (A·m²)

            % ================== MAGNET POSITION DATA ======================
            % Extract position and rotation data for all magnets

            num_magnets = size(epms.positions, 1);
            positions = epms.positions(:,1:2);  % Extract XY positions
            psai_vals = epms.psai;              % Rotation parameters

            % Calculate constant coefficient for all magnets:
            % C_i = (3μ₀/4π) * m_source * m_robot (precompute for each magnet)
            C_base = (3 * mu0) / (4 * pi) * m_source * m_robot;

            
            
            
            
            % ============== GENERATE FORCE FIELD FUNCTION =================
            % Create MATLAB function file with unrolled calculations

            fid = fopen('force_field_fast.m', 'w');

            % Function Header
            fprintf(fid, 'function [F_total, H] = force_field_fast(x, y, psai)\n');
            fprintf(fid, '%%FORCE_FIELD_FAST  Compute net magnetic force and Hessian\n');
            fprintf(fid, '%% [F_total, H] = force_field_fast(x, y, psai)\n');
            fprintf(fid, '%% Inputs:\n');
            fprintf(fid, '%%   x,y    - Robot position (meters)\n');
            fprintf(fid, '%%   psai   - Magnet rotation vector (radians)\n');
            fprintf(fid, '%% Outputs:\n');
            fprintf(fid, '%%   F_total - Net force [Fx, Fy, 0] (N)\n');
            fprintf(fid, '%%   H       - Hessian matrix [[Hxx, Hxy]; [Hxy, Hyy]] (N/m)\n\n');

            % Constants Section
            fprintf(fid, '%% ========== PHYSICAL CONSTANTS ==========\n');
            fprintf(fid, 'mu0 = %.12e;  %% Vacuum permeability (T·m/A)\n', mu0);
            fprintf(fid, 'm_robot = %.12e;  %% Robot dipole moment (A·m²)\n\n', m_robot);

            % Force/Hessian Initialization
            fprintf(fid, '%% ========== INITIALIZE ACCUMULATORS ==========\n');
            fprintf(fid, 'Fx_total = 0;  %% X-force accumulator\n');
            fprintf(fid, 'Fy_total = 0;  %% Y-force accumulator\n');
            fprintf(fid, 'Hxx_total = 0; %% Hessian XX component\n');
            fprintf(fid, 'Hxy_total = 0; %% Hessian XY component\n');
            fprintf(fid, 'Hyy_total = 0; %% Hessian YY component\n\n');

            % Individual Magnet Contributions
            for i = 1:num_magnets
                % Current magnet parameters
                px = positions(i,1);
                py = positions(i,2);
                C_i = C_base * cos(psai_vals(i));  % Precompute rotation component

                % Distance Calculations
                fprintf(fid, '%% ========== MAGNET %d @ [%.4f, %.4f] ==========\n', i, px, py);

                % Relative position components
                fprintf(fid, 'dx%d = x - %.12e;\n', i, px);
                fprintf(fid, 'dy%d = y - %.12e;\n', i, py);

                % Distance calculation
                fprintf(fid, 'r%d_sq = dx%d^2 + dy%d^2;\n', i, i, i);
                fprintf(fid, 'r%d = sqrt(r%d_sq);\n', i, i);
%                 fprintf(fid, 'if r%d == 0, continue; end  %% Skip singularity\n\n', i);

                % Force Components
                fprintf(fid, '%% ----- Force Contribution -----\n');
                fprintf(fid, 'r%d_5 = r%d^5;\n', i, i);  % r^5 term
                fprintf(fid, 'F_coeff%d = %.12e / r%d_5;\n', i, C_i, i);

                % X-component
                fprintf(fid, 'Fx_total = Fx_total + F_coeff%d * dx%d;\n', i, i);
                % Y-component
                fprintf(fid, 'Fy_total = Fy_total + F_coeff%d * dy%d;\n\n', i, i);

                % Hessian Components
                fprintf(fid, '%% ----- Hessian Contribution -----\n');
                fprintf(fid, 'r%d_7 = r%d_sq * r%d_5;  %% r^7 equivalent\n', i, i, i);

                % Common coefficient for Hessian terms
                fprintf(fid, 'H_coeff%d = %.12e / r%d_7;\n', i, -5*C_i, i);

                % XX component: (1/r^5 - 5*dx²/r^7)
                fprintf(fid, 'Hxx_total = Hxx_total + (%.12e/r%d_5) + H_coeff%d*dx%d^2;\n', C_i, i, i, i);

                % XY component: -5*dx*dy/r^7
                fprintf(fid, 'Hxy_total = Hxy_total + H_coeff%d*dx%d*dy%d;\n', i, i, i);

                % YY component: (1/r^5 - 5*dy²/r^7)
                fprintf(fid, 'Hyy_total = Hyy_total + (%.12e/r%d_5) + H_coeff%d*dy%d^2;\n\n', C_i, i, i, i);
            end

            % Final Output Assembly
            fprintf(fid, '%% ========== FINAL OUTPUT ASSEMBLY ==========\n');
            fprintf(fid, 'F_total = [Fx_total, Fy_total, 0];\n');
            fprintf(fid, 'H = [Hxx_total, Hxy_total;\n');
            fprintf(fid, '     Hxy_total, Hyy_total];\n');
            fprintf(fid, 'end\n');

            % Cleanup
            fclose(fid);
            fprintf('Optimized force field function with Hessian generated successfully.\n');
            
            
            
            
            
            
            % ============== GENERATE FORCE & GRADIENT FUNCTION ==============
            % Create MATLAB function file for force + gradients
            fid_grad = fopen('force_gradients_fast.m', 'w');

            % Function Header for Gradients
            fprintf(fid_grad, 'function [F_total, Grad] = force_gradients_fast(x, y, psai)\n');
            fprintf(fid_grad, '%%FORCE_GRADIENTS_FAST  Compute force and gradients w.r.t psai\n');
            fprintf(fid_grad, '%% [F_total, Grad] = force_gradients_fast(x, y, psai)\n');
            fprintf(fid_grad, '%% Grad(i) = [dFx/dpsai_i, dFy/dpsai_i]\n\n');

            % Constants Section (same as before)
            fprintf(fid_grad, 'mu0 = %.12e;\n', mu0);
            fprintf(fid_grad, 'm_robot = %.12e;\n\n', m_robot);

            % Initialize Accumulators
            fprintf(fid_grad, 'Fx = 0; Fy = 0;\n');
            fprintf(fid_grad, 'Grad = zeros(%d, 2); %% [dFx/dpsai_i; dFy/dpsai_i]\n\n', num_magnets);

            % Individual Magnet Contributions (Modified)
            for i = 1:num_magnets
                px = positions(i,1);
                py = positions(i,2);
                C_i = (3 * mu0 / (4 * pi)) * m_source * m_robot; % No cos(psai) here

                fprintf(fid_grad, '%% === Magnet %d @ [%.4f, %.4f] ===\n', i, px, py);

                % Relative position calculations
                fprintf(fid_grad, 'dx%d = x - %.12e;\n', i, px);
                fprintf(fid_grad, 'dy%d = y - %.12e;\n', i, py);
                fprintf(fid_grad, 'r%d = sqrt(dx%d^2 + dy%d^2);\n', i, i, i);
                %fprintf(fid_grad, 'if r%d == 0, continue; end\n\n', i);

                % Force components with dynamic psai
                fprintf(fid_grad, 'F_coeff%d = %.12e / r%d^5;\n', i, C_i, i);
                fprintf(fid_grad, 'Fx = Fx + F_coeff%d * dx%d * cos(psai(%d));\n', i, i, i);
                fprintf(fid_grad, 'Fy = Fy + F_coeff%d * dy%d * cos(psai(%d));\n\n', i, i, i);

                % Gradient components (analytical derivatives)
                fprintf(fid_grad, 'Grad(%d,1) = -F_coeff%d * dx%d * sin(psai(%d));\n', i, i, i, i);
                fprintf(fid_grad, 'Grad(%d,2) = -F_coeff%d * dy%d * sin(psai(%d));\n\n', i, i, i, i);
            end

            % Final Output Assembly
            fprintf(fid_grad, 'F_total = [Fx, Fy, 0];\n');
            fprintf(fid_grad, 'end\n');
            fclose(fid_grad);

            fprintf('Force + gradients function generated successfully.\n');
        end

        
        
        function saveObject(obj, filename)
            %SAVEOBJECT Save the current object to a MAT file.
            s = obj;
            save(filename, 's');
        end
    end
    

end
