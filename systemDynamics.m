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
m_mr = inputs.m_mr;
m_fp = inputs.m_fp;
i_fp = inputs.i_fp;
Cd = inputs.drag_coeff;
Mu = inputs.viscosity;
mu_0 = inputs.args.mu_0;
m_ = inputs.args.mr.m / 5e1; %%% note
threshold = inputs.threshold;
sigma = inputs.sigma;
epsilun = inputs.epsilun;
%
LJ_potential = @(r)4*epsilun*( (sigma/r)^12 - (sigma/r)^6 );
LJ_force     = @(r)4*epsilun*( -12*(sigma/r)^12/r + 6*(sigma/r)^6/r );
% LJ_force     = @(r)4*epsilun*( + 6*(sigma/r)^6/r );
% LJ_force     = @(r)0.5e-8;
%
D_e = epsilun;
r_e = sigma;
k_e = 100;
a = k_e/2/D_e;
MorsePotential = @(r) D_e * ( 1-exp(-a*(r-r_e)) )^2;
%
% sigma2 = 1 * 1e-5;
% epsilun2 = 0.01;
% LJ_potential2 = @(r)4*epsilun2*( (sigma2/r)^12 - (sigma2/r)^6 );
% LJ_force2     = @(r,epsilun)4*epsilun*( -12*(sigma2/r)^12/r + 6*(sigma2/r)^6/r );

% magnethic forces
for i=1:n
    F = force_field_symbolic(x_mr(i), y_mr(i), Psai);
    fx(i) = F(1);
    fy(i) = F(2);
end
% drag forces of MRs
for i=1:n
    fx(i) = fx(i) - Mu * inputs.k_mr(1,i) * xd_mr(i);
    fy(i) = fy(i) - Mu * inputs.k_mr(2,i) * yd_mr(i);
end
% drag forces of FPs
for cnt = 1:m
    fx(cnt+n) = fx(cnt+n) - Mu * inputs.k_mr(1,i) * xd_fp(cnt) * 5;
    fy(cnt+n) = fy(cnt+n) - Mu * inputs.k_mr(2,i) * yd_fp(cnt) * 5;
end
% interaction force between MRs
for i=1:n
    for j=1:n
        if i ~= j
            vector =  [x_mr(i) y_mr(i)] - [x_mr(j) y_mr(j)]; % direction from j to i
            mr_radius1 = inputs.r_mr(i);
            mr_radius2 = inputs.r_mr(j);
            r1 = norm(vector);
            r2 = norm(vector) - mr_radius1 - mr_radius2;
            omega = 3*mu_0*m_*m_/4/pi;
            force = omega / r1^4;
            normal_i = vector ./ norm(vector);
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
        mr_radius = inputs.r_mr(i);
        if fp.type == 1
            fp_radius = fp.radius;
            vector =  [x_fp(cnt) y_fp(cnt)] - [x_mr(i) y_mr(i)]; % direction from mr to fp
            r =  norm(vector) - fp_radius - mr_radius ;
            if r < threshold
                normal_i = vector ./ norm(vector);
                reactionForceMax = 5.001 * abs( dot([fx(i) fy(i)], normal_i) );
                force = min(LJ_force(r), reactionForceMax); %% note
                if force < 0
                    force = reactionForceMax;
                end
                fx_buff = fx_buff + normal_i(1)*force;
                fy_buff = fy_buff + normal_i(2)*force;
                %%%
                force;
                if force < 0
                    a=1;
                end
                force1x = fx(i);
                force1y = fy(i);
                forcex = fx_buff / force1x;
                forcey = fy_buff / force1y;
                %%%
            end
        elseif fp.type == 2
            vector =  [x_fp(cnt) y_fp(cnt)] - [x_mr(i) y_mr(i)];  % direction from mr to fp
            reactionForceMax = 10.001 * abs( dot([fx(i) fy(i)], normal_i) );
            fp_length = fp.length;
            theta_vec = atan2(vector(2), vector(1));
            theta_rec = t_fp(cnt);
            theta = theta_vec - theta_rec;
            r_x = abs(vector(1)) - fp_length/2 - mr_radius;
            r_y = abs(vector(2)) - fp_length/2 - mr_radius;
            r__ = norm(vector) - sqrt(2)/2*fp_length - mr_radius;
            r = max(r_x, r_y);
            if r__ < threshold
                force = min(LJ_force(r), reactionForceMax); %% note
                if force < 0
                    force = reactionForceMax;
                end
                if r_x < r_y
                    normal_i = [1 0] * sign(vector(1));
                else
                    normal_i = [0 1] * sign(vector(2));
                end
                fx_buff = fx_buff + normal_i(1)*force;
                fy_buff = fy_buff + normal_i(2)*force;
            end
        elseif fp.type == 3
            fp_points = fp.points + [x_fp(cnt); y_fp(cnt)];
            wall{1} = [fp_points(1:2,1); fp_points(1:2,end)];
            for  j=1:size(fp_points,2)-1
                wall{j+1} = [fp_points(1:2,j); fp_points(1:2,j+1)];
            end
            for  j=1:size(fp_points,2)
                [r,isInContact,normal_i] = calculateDistanceToWallLinear([x_mr(i); y_mr(i)], wall{j}, threshold + mr_radius);% direction wall to point
                r = r - mr_radius;
                if isInContact
                    force = min(LJ_force(r), contactForceaMax); %% note
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
%         if fx_buff > 0
%             fx_buff
%             fy_buff
%         end
    end
end
%
% interaction force between Walls and MRs
for wlt = 1:w
    for i=1:n
        mr_radius = inputs.r_mr(i);
        if inputs.walls(1,wlt) == 0 % wall type
            [r,isInContact,normal_i] = calculateDistanceToWallLinear([x_mr(i); y_mr(i)], inputs.walls(2:5,wlt), threshold + mr_radius);% direction wall to point
            r = r - mr_radius;
        elseif inputs.walls(1,wlt) == 1
            a=1;
        end
        if isInContact
            reactionForceMax = 10.001 * abs( dot([fx(i) fy(i)], normal_i) );
            force = min(LJ_force(r), reactionForceMax);
            if force < 0
                force = reactionForceMax;
            end
            fx(i) = fx(i) + normal_i(1)*force;
            fy(i) = fy(i) + normal_i(2)*force;
        end
    end
end
% interaction force between Walls and FPs   %% fp type ok shavad
% for wlt = 1:w
%     for cnt=1:m
%         fp = inputs.fps{cnt};
%         fp_radius = fp.radius;
%         if inputs.walls(1,wlt) == 0 % wall type
%             [r,isInContact,normal_i] = calculateDistanceToWallLinear([x_fp(cnt); y_fp(cnt)], inputs.walls(2:5,wlt), threshold + fp_radius);% direction wall to point
%             r = r - fp_radius;
%         elseif inputs.walls(1,wlt) == 1
%             a=1;
%         end
%         if isInContact
%             reactionForceMax = 10.001 * abs( dot([fx(cnt+n) fy(cnt+n)], normal_i) );
%             force = min(LJ_force(r), reactionForceMax);
%             if force < 0
%                 force = reactionForceMax;
%             end
%             fx(cnt+n) = fx(cnt+n) + normal_i(1)*force;
%             fy(cnt+n) = fy(cnt+n) + normal_i(2)*force;
%         end
%     end
% end



dydt(ind1) = y(ind3);
dydt(ind2) = y(ind4);
dydt(ind3) = fx(ind5,1)./m_mr';
dydt(ind4) = fy(ind5,1)./m_mr';
%
dydt(ind1_) = y(ind3_);
dydt(ind2_) = y(ind4_);
dydt(ind3_) = fx(ind5_,1)./m_fp';
dydt(ind4_) = fy(ind5_,1)./m_fp';
dydt(ind6_) = (mo(ind5_,1)-0.*Cd)./i_fp';
dydt(ind7_) = (mo(ind5_,1)-0.*Cd)./i_fp';



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