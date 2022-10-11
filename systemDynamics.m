function dydt = systemDynamics(t,y)
n = size(y,1)/4;
x_ = y(1:n,1);
y_ = y(n+1:2*n,1);
xd_ = y(2*n+1:3*n,1);
yd_ = y(3*n+1:4*n,1);
fx = zeros(n,1);
fy = zeros(n,1);
dydt = zeros(size(y,1),1);
%
Psai = psaiController(t);
mass = [
    4/3*pi*0.025^3*(3) / 1e3 %Rho_pla = 1.25 g/cm3 Rho_Neodymium = 7.5 g/cm3 r=0.025 cm
];
c = [
    16/3*0.01*(0.25/1e3)
];
mu_0 = 4*pi*1e-7;
M = 1.2706/mu_0;              % Magnetization   [A/m]    
L = 0.0004;
D = 0.0002;
m = M * (pi*D^2/4*L);
%
for i=1:n
    F = force_field_symbolic(x_(i), y_(i), Psai);
    fx(i) = F(1);
    fy(i) = F(2);
end
%
for i=1:n-1
    for j=i+1:n
        r = sqrt( (x_(i)-x_(j))^2 + (y_(i)-y_(j))^2 );
        omega = 3*mu_0*m*m/4/pi;
        normal_i = [x_(i)-x_(j) y_(i)-y_(j)];  % direction from j to i
        normal_i = normal_i ./ norm(normal_i);
        fx(i) = fx(i) + normal_i(1)*omega/r^4;
        fy(i) = fy(i) + normal_i(2)*omega/r^4;
        fx(j) = fx(j) - normal_i(1)*omega/r^4;
        fy(j) = fy(j) - normal_i(2)*omega/r^4;
    end
end

dydt(1:n) = y(2*n+1:3*n);
dydt(n+1:2*n) = y(3*n+1:4*n);
dydt(2*n+1:3*n) = (fx-y(2*n+1:3*n).*c)./mass;
dydt(3*n+1:4*n) = (fy-y(3*n+1:4*n).*c)./mass;