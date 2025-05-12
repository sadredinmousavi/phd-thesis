function [F_total, Grad] = force_gradients_fast(x, y, psai)
%FORCE_GRADIENTS_FAST  Compute force and gradients w.r.t psai
% [F_total, Grad] = force_gradients_fast(x, y, psai)
% Grad(i) = [dFx/dpsai_i, dFy/dpsai_i]

mu0 = 1.256637061436e-06;
m_robot = 1.250000000000e-04;

Fx = 0; Fy = 0;
Grad = zeros(8, 2); % [dFx/dpsai_i; dFy/dpsai_i]

% === Magnet 1 @ [0.2500, 0.0000] ===
dx1 = x - 2.500000000000e-01;
dy1 = y - 0.000000000000e+00;
r1 = sqrt(dx1^2 + dy1^2);
F_coeff1 = 3.103521390292e-10 / r1^5;
Fx = Fx + F_coeff1 * dx1 * cos(psai(1));
Fy = Fy + F_coeff1 * dy1 * cos(psai(1));

Grad(1,1) = -F_coeff1 * dx1 * sin(psai(1));
Grad(1,2) = -F_coeff1 * dy1 * sin(psai(1));

% === Magnet 2 @ [0.1768, 0.1768] ===
dx2 = x - 1.767766952966e-01;
dy2 = y - 1.767766952966e-01;
r2 = sqrt(dx2^2 + dy2^2);
F_coeff2 = 3.103521390292e-10 / r2^5;
Fx = Fx + F_coeff2 * dx2 * cos(psai(2));
Fy = Fy + F_coeff2 * dy2 * cos(psai(2));

Grad(2,1) = -F_coeff2 * dx2 * sin(psai(2));
Grad(2,2) = -F_coeff2 * dy2 * sin(psai(2));

% === Magnet 3 @ [0.0000, 0.2500] ===
dx3 = x - 1.530808498934e-17;
dy3 = y - 2.500000000000e-01;
r3 = sqrt(dx3^2 + dy3^2);
F_coeff3 = 3.103521390292e-10 / r3^5;
Fx = Fx + F_coeff3 * dx3 * cos(psai(3));
Fy = Fy + F_coeff3 * dy3 * cos(psai(3));

Grad(3,1) = -F_coeff3 * dx3 * sin(psai(3));
Grad(3,2) = -F_coeff3 * dy3 * sin(psai(3));

% === Magnet 4 @ [-0.1768, 0.1768] ===
dx4 = x - -1.767766952966e-01;
dy4 = y - 1.767766952966e-01;
r4 = sqrt(dx4^2 + dy4^2);
F_coeff4 = 3.103521390292e-10 / r4^5;
Fx = Fx + F_coeff4 * dx4 * cos(psai(4));
Fy = Fy + F_coeff4 * dy4 * cos(psai(4));

Grad(4,1) = -F_coeff4 * dx4 * sin(psai(4));
Grad(4,2) = -F_coeff4 * dy4 * sin(psai(4));

% === Magnet 5 @ [-0.2500, 0.0000] ===
dx5 = x - -2.500000000000e-01;
dy5 = y - 3.061616997868e-17;
r5 = sqrt(dx5^2 + dy5^2);
F_coeff5 = 3.103521390292e-10 / r5^5;
Fx = Fx + F_coeff5 * dx5 * cos(psai(5));
Fy = Fy + F_coeff5 * dy5 * cos(psai(5));

Grad(5,1) = -F_coeff5 * dx5 * sin(psai(5));
Grad(5,2) = -F_coeff5 * dy5 * sin(psai(5));

% === Magnet 6 @ [-0.1768, -0.1768] ===
dx6 = x - -1.767766952966e-01;
dy6 = y - -1.767766952966e-01;
r6 = sqrt(dx6^2 + dy6^2);
F_coeff6 = 3.103521390292e-10 / r6^5;
Fx = Fx + F_coeff6 * dx6 * cos(psai(6));
Fy = Fy + F_coeff6 * dy6 * cos(psai(6));

Grad(6,1) = -F_coeff6 * dx6 * sin(psai(6));
Grad(6,2) = -F_coeff6 * dy6 * sin(psai(6));

% === Magnet 7 @ [-0.0000, -0.2500] ===
dx7 = x - -4.592425496803e-17;
dy7 = y - -2.500000000000e-01;
r7 = sqrt(dx7^2 + dy7^2);
F_coeff7 = 3.103521390292e-10 / r7^5;
Fx = Fx + F_coeff7 * dx7 * cos(psai(7));
Fy = Fy + F_coeff7 * dy7 * cos(psai(7));

Grad(7,1) = -F_coeff7 * dx7 * sin(psai(7));
Grad(7,2) = -F_coeff7 * dy7 * sin(psai(7));

% === Magnet 8 @ [0.1768, -0.1768] ===
dx8 = x - 1.767766952966e-01;
dy8 = y - -1.767766952966e-01;
r8 = sqrt(dx8^2 + dy8^2);
F_coeff8 = 3.103521390292e-10 / r8^5;
Fx = Fx + F_coeff8 * dx8 * cos(psai(8));
Fy = Fy + F_coeff8 * dy8 * cos(psai(8));

Grad(8,1) = -F_coeff8 * dx8 * sin(psai(8));
Grad(8,2) = -F_coeff8 * dy8 * sin(psai(8));

F_total = [Fx, Fy, 0];
end
