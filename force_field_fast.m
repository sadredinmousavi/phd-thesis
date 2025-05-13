function [F_total, H] = force_field_fast(x, y, psai)
%FORCE_FIELD_FAST  Compute net magnetic force and Hessian
% [F_total, H] = force_field_fast(x, y, psai)
% Inputs:
%   x,y    - Robot position (meters)
%   psai   - Magnet rotation vector (radians)
% Outputs:
%   F_total - Net force [Fx, Fy, 0] (N)
%   H       - Hessian matrix [[Hxx, Hxy]; [Hxy, Hyy]] (N/m)

% ========== PHYSICAL CONSTANTS ==========
mu0 = 1.256637061436e-06;  % Vacuum permeability (T·m/A)
m_robot = 1.250000000000e-04;  % Robot dipole moment (A·m²)

% ========== INITIALIZE ACCUMULATORS ==========
Fx_total = 0;  % X-force accumulator
Fy_total = 0;  % Y-force accumulator
Hxx_total = 0; % Hessian XX component
Hxy_total = 0; % Hessian XY component
Hyy_total = 0; % Hessian YY component

% ========== MAGNET 1 @ [0.2500, 0.0000] ==========
dx1 = x - 2.500000000000e-01;
dy1 = y - 0.000000000000e+00;
r1_sq = dx1^2 + dy1^2;
r1 = sqrt(r1_sq);
% ----- Force Contribution -----
r1_5 = r1^5;
F_coeff1 = 3.103521390292e-10 / r1_5;
Fx_total = Fx_total + F_coeff1 * dx1 * cos(psai(1));
Fy_total = Fy_total + F_coeff1 * dy1 * cos(psai(1));

% ----- Hessian Contribution -----
r1_7 = r1_sq * r1_5;  % r^7 equivalent
H_coeff1 = -1.551760695146e-09 / r1_7;
Hxx_total = Hxx_total + ( (3.103521390292e-10/r1_5) + H_coeff1*dx1^2 ) * cos(psai(1));
Hxy_total = Hxy_total + ( H_coeff1*dx1*dy1 ) * cos(psai(1));
Hyy_total = Hyy_total + ( (3.103521390292e-10/r1_5) + H_coeff1*dy1^2 ) * cos(psai(1));

% ========== MAGNET 2 @ [0.1768, 0.1768] ==========
dx2 = x - 1.767766952966e-01;
dy2 = y - 1.767766952966e-01;
r2_sq = dx2^2 + dy2^2;
r2 = sqrt(r2_sq);
% ----- Force Contribution -----
r2_5 = r2^5;
F_coeff2 = 3.103521390292e-10 / r2_5;
Fx_total = Fx_total + F_coeff2 * dx2 * cos(psai(2));
Fy_total = Fy_total + F_coeff2 * dy2 * cos(psai(2));

% ----- Hessian Contribution -----
r2_7 = r2_sq * r2_5;  % r^7 equivalent
H_coeff2 = -1.551760695146e-09 / r2_7;
Hxx_total = Hxx_total + ( (3.103521390292e-10/r2_5) + H_coeff2*dx2^2 ) * cos(psai(2));
Hxy_total = Hxy_total + ( H_coeff2*dx2*dy2 ) * cos(psai(2));
Hyy_total = Hyy_total + ( (3.103521390292e-10/r2_5) + H_coeff2*dy2^2 ) * cos(psai(2));

% ========== MAGNET 3 @ [0.0000, 0.2500] ==========
dx3 = x - 1.530808498934e-17;
dy3 = y - 2.500000000000e-01;
r3_sq = dx3^2 + dy3^2;
r3 = sqrt(r3_sq);
% ----- Force Contribution -----
r3_5 = r3^5;
F_coeff3 = 3.103521390292e-10 / r3_5;
Fx_total = Fx_total + F_coeff3 * dx3 * cos(psai(3));
Fy_total = Fy_total + F_coeff3 * dy3 * cos(psai(3));

% ----- Hessian Contribution -----
r3_7 = r3_sq * r3_5;  % r^7 equivalent
H_coeff3 = -1.551760695146e-09 / r3_7;
Hxx_total = Hxx_total + ( (3.103521390292e-10/r3_5) + H_coeff3*dx3^2 ) * cos(psai(3));
Hxy_total = Hxy_total + ( H_coeff3*dx3*dy3 ) * cos(psai(3));
Hyy_total = Hyy_total + ( (3.103521390292e-10/r3_5) + H_coeff3*dy3^2 ) * cos(psai(3));

% ========== MAGNET 4 @ [-0.1768, 0.1768] ==========
dx4 = x - -1.767766952966e-01;
dy4 = y - 1.767766952966e-01;
r4_sq = dx4^2 + dy4^2;
r4 = sqrt(r4_sq);
% ----- Force Contribution -----
r4_5 = r4^5;
F_coeff4 = 3.103521390292e-10 / r4_5;
Fx_total = Fx_total + F_coeff4 * dx4 * cos(psai(4));
Fy_total = Fy_total + F_coeff4 * dy4 * cos(psai(4));

% ----- Hessian Contribution -----
r4_7 = r4_sq * r4_5;  % r^7 equivalent
H_coeff4 = -1.551760695146e-09 / r4_7;
Hxx_total = Hxx_total + ( (3.103521390292e-10/r4_5) + H_coeff4*dx4^2 ) * cos(psai(4));
Hxy_total = Hxy_total + ( H_coeff4*dx4*dy4 ) * cos(psai(4));
Hyy_total = Hyy_total + ( (3.103521390292e-10/r4_5) + H_coeff4*dy4^2 ) * cos(psai(4));

% ========== MAGNET 5 @ [-0.2500, 0.0000] ==========
dx5 = x - -2.500000000000e-01;
dy5 = y - 3.061616997868e-17;
r5_sq = dx5^2 + dy5^2;
r5 = sqrt(r5_sq);
% ----- Force Contribution -----
r5_5 = r5^5;
F_coeff5 = 3.103521390292e-10 / r5_5;
Fx_total = Fx_total + F_coeff5 * dx5 * cos(psai(5));
Fy_total = Fy_total + F_coeff5 * dy5 * cos(psai(5));

% ----- Hessian Contribution -----
r5_7 = r5_sq * r5_5;  % r^7 equivalent
H_coeff5 = -1.551760695146e-09 / r5_7;
Hxx_total = Hxx_total + ( (3.103521390292e-10/r5_5) + H_coeff5*dx5^2 ) * cos(psai(5));
Hxy_total = Hxy_total + ( H_coeff5*dx5*dy5 ) * cos(psai(5));
Hyy_total = Hyy_total + ( (3.103521390292e-10/r5_5) + H_coeff5*dy5^2 ) * cos(psai(5));

% ========== MAGNET 6 @ [-0.1768, -0.1768] ==========
dx6 = x - -1.767766952966e-01;
dy6 = y - -1.767766952966e-01;
r6_sq = dx6^2 + dy6^2;
r6 = sqrt(r6_sq);
% ----- Force Contribution -----
r6_5 = r6^5;
F_coeff6 = 3.103521390292e-10 / r6_5;
Fx_total = Fx_total + F_coeff6 * dx6 * cos(psai(6));
Fy_total = Fy_total + F_coeff6 * dy6 * cos(psai(6));

% ----- Hessian Contribution -----
r6_7 = r6_sq * r6_5;  % r^7 equivalent
H_coeff6 = -1.551760695146e-09 / r6_7;
Hxx_total = Hxx_total + ( (3.103521390292e-10/r6_5) + H_coeff6*dx6^2 ) * cos(psai(6));
Hxy_total = Hxy_total + ( H_coeff6*dx6*dy6 ) * cos(psai(6));
Hyy_total = Hyy_total + ( (3.103521390292e-10/r6_5) + H_coeff6*dy6^2 ) * cos(psai(6));

% ========== MAGNET 7 @ [-0.0000, -0.2500] ==========
dx7 = x - -4.592425496803e-17;
dy7 = y - -2.500000000000e-01;
r7_sq = dx7^2 + dy7^2;
r7 = sqrt(r7_sq);
% ----- Force Contribution -----
r7_5 = r7^5;
F_coeff7 = 3.103521390292e-10 / r7_5;
Fx_total = Fx_total + F_coeff7 * dx7 * cos(psai(7));
Fy_total = Fy_total + F_coeff7 * dy7 * cos(psai(7));

% ----- Hessian Contribution -----
r7_7 = r7_sq * r7_5;  % r^7 equivalent
H_coeff7 = -1.551760695146e-09 / r7_7;
Hxx_total = Hxx_total + ( (3.103521390292e-10/r7_5) + H_coeff7*dx7^2 ) * cos(psai(7));
Hxy_total = Hxy_total + ( H_coeff7*dx7*dy7 ) * cos(psai(7));
Hyy_total = Hyy_total + ( (3.103521390292e-10/r7_5) + H_coeff7*dy7^2 ) * cos(psai(7));

% ========== MAGNET 8 @ [0.1768, -0.1768] ==========
dx8 = x - 1.767766952966e-01;
dy8 = y - -1.767766952966e-01;
r8_sq = dx8^2 + dy8^2;
r8 = sqrt(r8_sq);
% ----- Force Contribution -----
r8_5 = r8^5;
F_coeff8 = 3.103521390292e-10 / r8_5;
Fx_total = Fx_total + F_coeff8 * dx8 * cos(psai(8));
Fy_total = Fy_total + F_coeff8 * dy8 * cos(psai(8));

% ----- Hessian Contribution -----
r8_7 = r8_sq * r8_5;  % r^7 equivalent
H_coeff8 = -1.551760695146e-09 / r8_7;
Hxx_total = Hxx_total + ( (3.103521390292e-10/r8_5) + H_coeff8*dx8^2 ) * cos(psai(8));
Hxy_total = Hxy_total + ( H_coeff8*dx8*dy8 ) * cos(psai(8));
Hyy_total = Hyy_total + ( (3.103521390292e-10/r8_5) + H_coeff8*dy8^2 ) * cos(psai(8));

% ========== FINAL OUTPUT ASSEMBLY ==========
F_total = [Fx_total, Fy_total, 0];
H = [Hxx_total, Hxy_total;
     Hxy_total, Hyy_total];
end
