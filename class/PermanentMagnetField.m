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
                end
                obj.magnets = magnetsStruct;
            else
                obj.magnets = [];
            end
        end
        
        function [Bx, By, Bz] = calculateField(obj, x, y, z)
            %CALCULATEFIELD Calculate the magnetic field at point(s) (x,y,z)
            % using the dipole approximation.
            % x, y, and z must be arrays of the same size.
            Bx = zeros(size(x));
            By = zeros(size(y));
            Bz = zeros(size(z));
            for i = 1:numel(obj.magnets)
                mag = obj.magnets(i);
                [dBx, dBy, dBz] = obj.dipoleField(x, y, z, mag);
                Bx = Bx + dBx;
                By = By + dBy;
                Bz = Bz + dBz;
            end
        end
        
        
        function F_total = calculateForceOnMicrorobot(obj, P_robot, m_robot)
            % CALCULATEFORCEONMICROROBOT computes the net force on a microrobot
            % (with dipole moment m_robot located at P_robot) due to all the
            % permanent magnets (stored in obj.magnets).
            %
            % Inputs:
            %   P_robot : 1x3 position vector of the microrobot.
            %   m_robot : 1x3 dipole moment vector of the microrobot.
            %
            % Output:
            %   F_total : 1x3 force vector [N] acting on the microrobot.

            F_total = [0, 0, 0];  % Initialize net force
            for i = 1:numel(obj.magnets)
                mag = obj.magnets(i);
                % For this source magnet, use its stored position.
                P_source = mag.position;

                % Compute effective dipole moment for the source magnet.
                % If dimensions are provided, factor them in.
                if isfield(mag, 'dimensions') && ~isempty(mag.dimensions)
                    m_source = mag.Br / obj.mu_0 * prod(mag.dimensions) * mag.orientation;
                else
                    m_source = mag.Br / obj.mu_0 * mag.orientation;
                end

                % Calculate force from this magnet using our static method.
                f_i = PermanentMagnetField.calculateDipoleForce(obj.mu_0, P_source, m_source, P_robot, m_robot);
                F_total = F_total + f_i;
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
        
        
        function makeForceFunctionHradCode(obj)
            % makeForceFunctionFast Dynamically generates a MATLAB function file 
            % named "force_field_fast.m" that computes the in-plane force on a 
            % microrobot.
            %
            % The method reads the permanent magnets (epms) and the robot magnet 
            % parameters (pm_robot) from obj.defs. Each permanent magnet is assumed 
            % to have the fields:
            %      position    -> 1x3 vector [x,y,z] (only x and y are used)
            %      Br          -> remnant flux density [A/m]
            %      dimensions  -> 1x3 vector [a,b,c] in meters
            %      orientation -> 1x3 vector (not used in the in-plane force)
            %      psai        -> (optional) fixed rotation parameter (radians, default = 0)
            %
            % The pm_robot structure should include:
            %      M           -> robot magnetization [A/m]
            %      dimensions  -> 1x3 vector representing the magnet volume [a,b,c]
            %
            % This function precomputes constants and then “unrolls” the loop over
            % all magnets so that the generated function runs without loops or
            % symbolic evaluations.

            %%% Read parameters from defs
            % (Adjust these field names if your definitions structure uses different names.)
            epms = obj.defs.epms;
            pm_robot = obj.defs.microRobots.defaultMagnetParams;

            % Compute the robot's scalar dipole moment (m_robot)
            V_robot = prod(pm_robot.dimensions);     % volume of robot magnet [m^3]
            m_robot = pm_robot.M * V_robot;            % scalar robot moment

            % Number of magnets:
            Nmag = size(epms.positions, 1);

            % Preallocate arrays to hold per-magnet data.
            C_arr   = zeros(Nmag, 1);  % Coefficient for each magnet
            pos_arr = zeros(Nmag, 2);  % XY positions for each magnet
            psi_arr = zeros(Nmag, 1);  % Fixed magnet rotation parameters

            % The magnet volume and Br are assumed to be stored in epms.inputs.
            % (They are assumed uniform for all magnets here.)
            V_source = epms.inputs.V;
            Br = epms.inputs.Br;
            m_source = (Br / obj.mu_0) * V_source;  % Scalar dipole moment for each magnet

            for i = 1:Nmag
                pos_arr(i, :) = epms.positions(i, 1:2);  % Use only x and y
                psi_arr(i) = epms.psai(i);
                % Compute constant coefficient for magnet i:
                % C_i = (3*mu0/(4*pi)) * (m_source * m_robot)
                C_arr(i) = (3 * obj.mu_0 / (4*pi)) * (m_source * m_robot);
            end

            %%% Generate the file "force_field_fast.m"
            fid = fopen('force_field_fast.m', 'wt');
            if fid == -1
                error('Could not open file for writing.');
            end

            % Write file header and function signature.
            fprintf(fid, 'function F_total = force_field_fast(x,y,psai)\n');
            fprintf(fid, '%% force_field_fast  Computes the net in-plane force on a microrobot\n');
            fprintf(fid, '%% at position (x,y) with dipole orientation psai (radians).\n');
            fprintf(fid, '%% Auto-generated by Definitions.makeForceFunctionFast\n\n');

            % Write the constant mu0 and the robot dipole moment.
            fprintf(fid, 'mu0 = %g; %% Read from obj.mu_0\n', obj.mu_0);
            fprintf(fid, '%% Robot scalar dipole moment (m_robot):\n');
            fprintf(fid, 'm_robot = %g; %% (Computed from pm_robot.M and pm_robot.dimensions)\n\n', m_robot);

            % Initialize net force accumulators.
            fprintf(fid, '%% Initialize force accumulators (Fz = 0 since force is in-plane)\n');
            fprintf(fid, 'Fx_total = 0;\n');
            fprintf(fid, 'Fy_total = 0;\n\n');

            % Unroll the loop: For each magnet, write its contribution.
            for i = 1:Nmag
                p_x = pos_arr(i,1);
                p_y = pos_arr(i,2);
                C_i = C_arr(i);

                fprintf(fid, '%% Magnet %d (position: [%g, %g])\n', i, p_x, p_y);
                fprintf(fid, 'r%d = sqrt((x - (%.10e))^2 + (y - (%g))^2);\n', i, p_x, p_y);
                fprintf(fid, 'r%d_hat_x = (x - (%.10e)) / r%d;\n', i, p_x, i);
                fprintf(fid, 'r%d_hat_y = (y - (%.10e)) / r%d;\n', i, p_y, i);
                fprintf(fid, 'F%d_x = (%.10e / (r%d^4)) * cos(psai(%d)) * r%d_hat_x;\n', i, C_i, i, i, i);
                fprintf(fid, 'F%d_y = (%.10e / (r%d^4)) * cos(psai(%d)) * r%d_hat_y;\n', i, C_i, i, i, i);
                fprintf(fid, 'Fx_total = Fx_total + F%d_x;\n', i);
                fprintf(fid, 'Fy_total = Fy_total + F%d_y;\n\n', i);
            end

            % Sum contributions and form the output force vector.
            fprintf(fid, 'F_total = [Fx_total, Fy_total, 0];\n');
            fprintf(fid, 'end\n');

            fclose(fid);
            fprintf('force_field_fast.m has been generated successfully.\n');
        end

        
        
        function saveObject(obj, filename)
            %SAVEOBJECT Save the current object to a MAT file.
            s = obj;
            save(filename, 's');
        end
    end
    
    
    
    
    
    
    methods (Static)
        function obj = loadObject(filename)
            %LOADOBJECT Load a PermanentMagnetField object from a MAT file.
            tmp = load(filename, 's');
            obj = tmp.s;
        end
        
        function B = compiledDipoleCalc(r, mu_0, magPositions, baseDipoleMoments, psiVec)
            %COMPILEDDIPOLECALC Fast dipole field calculation.
            %   r - 1×3 observation point.
            %   psiVec - an N×1 vector of ψ values (radians) for each magnet.
            %   magPositions - N×3 array of magnet positions.
            %   baseDipoleMoments - N×3 array of base dipole moments.
            %
            % Computes each magnet's effective dipole moment by rotating
            % baseDipoleMoments(j,:) about the fixed z-axis.
            %
            % Returns B, a 3×1 magnetic field vector.
            B = [0, 0, 0];
            N = size(magPositions, 1);
            for j = 1:N
                dr = r - magPositions(j,:);
                r_norm = norm(dr);
                if r_norm == 0
                    continue;  % avoid singularity
                end
                r_hat = dr / r_norm;
                psi = psiVec(j);
                % Rotate about the z-axis using raw sine and cosine:
                Rz = [cos(psi), -sin(psi), 0;
                      sin(psi),  cos(psi), 0;
                      0,         0,        1];
                m_eff = Rz * baseDipoleMoments(j,:)';
                m_dot_r = dot(m_eff, r_hat);
                B_contrib = (mu_0/(4*pi)) / (r_norm^3) * (3*m_dot_r*r_hat - m_eff');
                B = B + B_contrib;
            end
            B = B(:);  % Return as a 3×1 column vector.
        end
        
        function f = calculateDipoleForce(mu_0, P_source, m_source, P_robot, m_robot)
            % CALCULATEDIPOLEFORCE calculates the force on a target dipole (m_robot at P_robot)
            % produced by a source dipole (m_source at P_source) using the dipole-dipole force law.
            %
            % Inputs:
            %   mu0      - Permeability of free space.
            %   P_source - 1x3 position vector of the source magnet.
            %   m_source - 1x3 dipole moment of the source magnet.
            %   P_robot  - 1x3 position vector of the microrobot (target dipole).
            %   m_robot  - 1x3 dipole moment of the microrobot.
            %
            % Output:
            %   f        - 1x3 force vector acting on the microrobot.

            r = P_robot - P_source;
            r_norm = norm(r);
            if r_norm == 0
                f = [0, 0, 0];
                return;
            end
            r_hat = r / r_norm;
            
            % Ensure r_hat, m_robot, and m_source are all row vectors
            r_hat = r_hat(:)';
            m_robot = m_robot(:)';
            m_source = m_source(:)';
            
            dot_r_mrobot = dot(r_hat, m_robot);
            dot_r_msource = dot(r_hat, m_source);
            
            term = (dot_r_mrobot) * m_source + (dot_r_msource) * m_robot - 5 * (dot_r_msource) * (dot_r_mrobot) * r_hat;
            f = (3 * mu_0 / (4 * pi)) * term / (r_norm^3);
        end
        
        
    
    end
    
    methods (Access = private)
        function [Bx, By, Bz] = dipoleField(obj, x, y, z, magnet)
            %DIPOLEFIELD Compute the dipole field for a single magnet.
            % Uses the dipole formula:
            %   B = (mu0/(4*pi*r^3)) * (3*(m·r_hat)*r_hat - m)
            %
            % where r = [x,y,z] - magnet.position.
            if isempty(magnet.position)
                error('Magnet position must be provided.');
            end
            m = magnet.Br / obj.mu_0 * magnet.orientation;
            Bx = zeros(size(x));
            By = zeros(size(y));
            Bz = zeros(size(z));
            for i = 1:numel(x)
                r_vec = [x(i)-magnet.position(1), ...
                         y(i)-magnet.position(2), ...
                         z(i)-magnet.position(3)];
                r_norm = norm(r_vec);
                if r_norm > 0
                    r_hat = r_vec / r_norm;
                    B_local = (obj.mu0/(4*pi*r_norm^3)) * (3*dot(m, r_hat)*r_hat - m);
                    Bx(i) = B_local(1);
                    By(i) = B_local(2);
                    Bz(i) = B_local(3);
                end
            end
        end
    end
end
