function dydt = systemDynamics(t,y, inputs)
n = inputs.mr_num; % n = size(y,1)/4;
m = inputs.fp_num;
w = size(inputs.walls,2);
%
x_mr = y(1:n,1);
y_mr = y(n+1:2*n,1);
xd_mr = y(2*n+1:3*n,1);
yd_mr = y(3*n+1:4*n,1);
%
x_fp = y(4*n+1:4*n+m,1);
y_fp = y(4*n+m+1:4*n+2*m,1);
xd_fp = y(4*n+2*m+1:4*n+3*m,1);
yd_fp = y(4*n+3*m+1:4*n+4*m,1);
%
fx = zeros(n+m,1);
fy = zeros(n+m,1);
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
m_ = M * (pi*D^2/4*L);
%
sigma = 1 * 1e-5;
epsilun = 1000;
LJ_potential = @(r)4*epsilun*( (sigma/r)^12 - (sigma/r)^6 );
LJ_force     = @(r)4*epsilun*( -12*(sigma/r)^12/r + 6*(sigma/r)^6/r );
%
for i=1:n
    F = force_field_symbolic(x_mr(i), y_mr(i), Psai);
    fx(i) = F(1);
    fy(i) = F(2);
end
% interaction force between MRs
for i=1:n
    for j=i:n
        if i ~= j
            r = sqrt( (x_mr(i)-x_mr(j))^2 + (y_mr(i)-y_mr(j))^2 );
            omega = 3*mu_0*m_*m_/4/pi;
            force = omega / r^4;
            normal_i = [x_mr(i)-x_mr(j) y_mr(i)-y_mr(j)];  % direction from j to i
            normal_i = normal_i ./ norm(normal_i);
            fx(i) = fx(i) + normal_i(1)*force;
            fy(i) = fy(i) + normal_i(2)*force;
            fx(j) = fx(j) - normal_i(1)*force;
            fy(j) = fy(j) - normal_i(2)*force;
        end
    end
end
% interaction force between FPs
for i=1:m
    for j=i:m
        if i ~= j
%             r = sqrt( (x_fp(i)-x_fp(j))^2 + (y_fp(i)-y_fp(j))^2 ) - r_fp(i) - r_fp(j);
            r = sqrt( (x_fp(i)-x_fp(j))^2 + (y_fp(i)-y_fp(j))^2 ) - 0.02 - 0.02;
            force = LJ_force(r);
            normal_i = [x_mr(i)-x_mr(j) y_mr(i)-y_mr(j)];  % direction from j to i
            normal_i = normal_i ./ norm(normal_i);
            fx(i) = fx(i) + normal_i(1)*force;
            fy(i) = fy(i) + normal_i(2)*force;
            fx(j) = fx(j) - normal_i(1)*force;
            fy(j) = fy(j) - normal_i(2)*force;
        end
    end
end
% interaction force between FPs and MRs
for cnt = 1:m
    for i=1:n
%         r = sqrt( (x_fp(cnt)-x_mr(i))^2 + (y_fp(cnt)-y_mr(i))^2 ) - r_fp(cnt) - r_mr(i);
        r = sqrt( (x_fp(cnt)-x_mr(i))^2 + (y_fp(cnt)-y_mr(i))^2 ) - 0.002 - 0.001;
        force = LJ_force(r);
        normal_i = [x_fp(cnt)-x_mr(i) y_fp(cnt)-y_mr(i)];  % direction from j to i
        normal_i = normal_i ./ norm(normal_i);
        fx(cnt+n) = fx(cnt+n) + normal_i(1)*force;
        fy(cnt+n) = fy(cnt+n) + normal_i(2)*force;
        fx(i) = fx(i) - normal_i(1)*force;
        fy(i) = fy(i) - normal_i(2)*force;
    end
end
%
% interaction force between Walls and MRs
for wlt = 1:w
    for i=1:n
        if inputs.walls(1,wlt) == 0
            threshold = 0.002;
            [r,isInContact,normal_i] = calculateDistanceToWallLinear([x_mr(i); y_mr(i)], inputs.walls(2:5,wlt), threshold);% direction wall to point
        elseif inputs.walls(1,wlt) == 1
            a=1;
        end
        if isInContact
            force = LJ_force(r);
            fx(i) = fx(i) + normal_i(1)*force;
            fy(i) = fy(i) + normal_i(2)*force;
        end
    end
end




% dydt(1:n) = y(2*n+1:3*n);
% dydt(n+1:2*n) = y(3*n+1:4*n);
% dydt(2*n+1:3*n) = (fx(1:n,1)-y(2*n+1:3*n).*c)./mass;
% dydt(3*n+1:4*n) = (fy(1:n,1)-y(3*n+1:4*n).*c)./mass;

ind1 = 1:n;
ind2 = n+1:2*n;
ind3 = 2*n+1:3*n;
ind4 = 3*n+1:4*n;
ind5 = 1:n;
%
dydt(ind1) = y(ind3);
dydt(ind2) = y(ind4);
dydt(ind3) = (fx(ind5,1)-y(ind3).*c)./mass;
dydt(ind4) = (fy(ind5,1)-y(ind4).*c)./mass;
%
ind1 = (1:m) + 4*n;
ind2 = (m+1:2*m) + 4*n;
ind3 = (2*m+1:3*m) + 4*n;
ind4 = (3*m+1:4*m) + 4*n;
ind5 = (1:m)+n;
%
dydt(ind1) = y(ind3);
dydt(ind2) = y(ind4);
dydt(ind3) = (fx(ind5,1)-y(ind3).*c)./mass;
dydt(ind4) = (fy(ind5,1)-y(ind4).*c)./mass;