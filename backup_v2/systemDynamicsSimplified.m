function dydt = systemDynamicsSimplified(t,y,Psai)
n = size(y,1)/4;
x_ = y(1:n,1);
y_ = y(n+1:2*n,1);
% xd_ = y(2*n+1:3*n,1);
% yd_ = y(3*n+1:4*n,1);
fx = zeros(n,1);
fy = zeros(n,1);
dydt = zeros(size(y,1),1);
%
% Psai = psaiController(t);
mass = [
    4/3*pi*0.025^3*(3) / 1e3 %Rho_pla = 1.25 g/cm3 Rho_Neodymium = 7.5 g/cm3 r=0.025 cm
];
c = [
    16/3*0.01*(0.25/1e3)
];
%
for i=1:n
    F = force_field_symbolic(x_(i), y_(i), Psai);
    fx(i) = F(1);
    fy(i) = F(2);
end
%

dydt(1:n) = y(2*n+1:3*n);
dydt(n+1:2*n) = y(3*n+1:4*n);
dydt(2*n+1:3*n) = (fx-y(2*n+1:3*n).*c)./mass;
dydt(3*n+1:4*n) = (fy-y(3*n+1:4*n).*c)./mass;