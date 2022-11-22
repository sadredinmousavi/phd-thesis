function dydt = systemDynamics(t,y, inputs)

n = inputs.mr_num; % n = size(y,1)/4;
m = inputs.fp_num;
w = size(inputs.walls,2);
[ind1, ind2, ind3, ind4, ind5] = makeIndiceForMRs(n);
[ind1_, ind2_, ind3_, ind4_, ind5_, ind6_, ind7_] = makeIndiceForFPs(n, m);
%
x_mr = y(ind1);
y_mr = y(ind2);
xd_mr = y(ind3);
yd_mr = y(ind4);
%
x_fp = y(ind1_);
y_fp = y(ind2_);
xd_fp = y(ind3_);
yd_fp = y(ind4_);
t_fp = y(ind6_); % stands for theta
td_fp = y(ind7_);
%
fx = zeros(n+m,1);
fy = zeros(n+m,1);
mo = zeros(n+m,1);
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
epsilun = 0.01;
LJ_potential = @(r)4*epsilun*( (sigma/r)^12 - (sigma/r)^6 );
LJ_force     = @(r)4*epsilun*( -12*(sigma/r)^12/r + 6*(sigma/r)^6/r );
fplot(LJ_potential,[sigma/10,3*sigma]);
fplot(LJ_force,[0.000001,0.0005])
%
sigma2 = 1 * 1e-5;
epsilun2 = 10;
LJ_potential2 = @(r)4*epsilun2*( (sigma2/r)^12 - (sigma2/r)^6 );
LJ_force2     = @(r)4*epsilun2*( -12*(sigma2/r)^12/r + 6*(sigma2/r)^6/r );
%
for i=1:n
    F = force_field_symbolic(x_mr(i), y_mr(i), Psai);
    fx(i) = F(1);
    fy(i) = F(2);
end
% interaction force between MRs
for i=1:n
    for j=1:n
        if i ~= j
            mr_radius_1 = inputs.r_mr(i);
            mr_radius_2 = inputs.r_mr(j);
            r = sqrt( (x_mr(i)-x_mr(j))^2 + (y_mr(i)-y_mr(j))^2 ) - mr_radius_1 - mr_radius_2;
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
    for j=1:m
        if i ~= j
%             r = sqrt( (x_fp(i)-x_fp(j))^2 + (y_fp(i)-y_fp(j))^2 ) - r_fp(i) - r_fp(j);
%             r = sqrt( (x_fp(i)-x_fp(j))^2 + (y_fp(i)-y_fp(j))^2 ) - 0.02 - 0.02;
%             force = LJ_force(r);
%             normal_i = [x_mr(i)-x_mr(j) y_mr(i)-y_mr(j)];  % direction from j to i
%             normal_i = normal_i ./ norm(normal_i);
%             fx(n+i) = fx(n+i) + normal_i(1)*force;
%             fy(n+i) = fy(n+i) + normal_i(2)*force;
%             fx(n+j) = fx(n+j) - normal_i(1)*force;
%             fy(n+j) = fy(n+j) - normal_i(2)*force;
        end
    end
end
% interaction force between FPs and MRs
for cnt = 1:m
    for i = 1:n
        fp = inputs.fps{cnt};
        fx_buff = 0;
        fy_buff = 0;
        mo_buff = 0;
        %
        mr_radius = inputs.r_mr(i);
        threshold = 3*mr_radius;
        if fp.type == 1
            fp_radius = fp.radius;
            r = sqrt( (x_fp(cnt)-x_mr(i))^2 + (y_fp(cnt)-y_mr(i))^2 ) - fp_radius - mr_radius;
            if r < threshold
                force = LJ_force(r);
                normal_i = [x_fp(cnt)-x_mr(i) y_fp(cnt)-y_mr(i)];  % direction from j to i
                normal_i = normal_i ./ norm(normal_i);
                fx_buff = fx_buff + normal_i(1)*force;
                fy_buff = fy_buff + normal_i(2)*force;
            end
        elseif fp.type == 2
            fp_points = fp.points + [x_fp(cnt); y_fp(cnt)];
            wall{1} = [fp_points(1:2,1); fp_points(1:2,end)];
            for  j=1:size(fp_points,2)-1
                wall{j+1} = [fp_points(1:2,j); fp_points(1:2,j+1)];
            end
            for  j=1:size(fp_points,2)
                [r,isInContact,normal_i] = calculateDistanceToWallLinear([x_mr(i); y_mr(i)], wall{j}, threshold);% direction wall to point
                r = r - mr_radius;
                if isInContact
                    force = LJ_force2(r);
                    fx_buff = fx_buff + normal_i(1)*force;
                    fy_buff = fy_buff + normal_i(2)*force;
%                     mo_buff = mo_buff + normal_i(2)*force;
                end
                if norm(wall{j})>1
                    a=1;
%                     t
                end
            end
    
        end
        fx(cnt+n) = fx(cnt+n) + fx_buff;
        fy(cnt+n) = fy(cnt+n) + fy_buff;
        fx(i) = fx(i) - fx_buff;
        fy(i) = fy(i) - fy_buff;
        if fx_buff > 0
            fx_buff
            fy_buff
        end
    end
end
%
% interaction force between Walls and MRs
for wlt = 1:w
    for i=1:n
        mr_radius = inputs.r_mr(i);
        if inputs.walls(1,wlt) == 0 % wall type
            threshold = 3*mr_radius;
            [r,isInContact,normal_i] = calculateDistanceToWallLinear([x_mr(i); y_mr(i)], inputs.walls(2:5,wlt), threshold);% direction wall to point
            r = r - mr_radius;
        elseif inputs.walls(1,wlt) == 1
            a=1;
        end
        if isInContact
            force = LJ_force(r);
            if force>0
                force
            end
            fx(i) = fx(i) + normal_i(1)*force;
            fy(i) = fy(i) + normal_i(2)*force;
        end
    end
end


dydt(ind1) = y(ind3);
dydt(ind2) = y(ind4);
dydt(ind3) = (fx(ind5,1)-y(ind3).*c)./mass;
dydt(ind4) = (fy(ind5,1)-y(ind4).*c)./mass;
%
dydt(ind1_) = y(ind3_);
dydt(ind2_) = y(ind4_);
dydt(ind3_) = (fx(ind5_,1)-y(ind3_).*c)./mass;
dydt(ind4_) = (fy(ind5_,1)-y(ind4_).*c)./mass;
dydt(ind6_) = (mo(ind5_,1)-0.*c)./mass;
dydt(ind7_) = (mo(ind5_,1)-0.*c)./mass;



function [ind1, ind2, ind3, ind4, ind5] = makeIndiceForMRs(n)
    ind1 = 1:n;
    ind2 = n+1:2*n;
    ind3 = 2*n+1:3*n;
    ind4 = 3*n+1:4*n;
    ind5 = 1:n;
end
%
function [ind1, ind2, ind3, ind4, ind5, ind6, ind7] = makeIndiceForFPs(n, m)
    ind1 = (1:m) + 4*n;
    ind2 = (m+1:2*m) + 4*n;
    ind3 = (2*m+1:3*m) + 4*n;
    ind4 = (3*m+1:4*m) + 4*n;
    ind5 = (1:m)+n;
    % note that 4*n must be changed if MRs dof changes
    ind6 = (4*m+1:5*m) + 4*n;
    ind7 = (5*m+1:6*m) + 4*n;
end
end