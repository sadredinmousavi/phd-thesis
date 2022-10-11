%% Cylinder magnetic field visualization
 % This script is a DEMO for the visualization of the magnetic flux density 
 % field lines of an axially magnetized cylinder. 

clc, close all, clearvars
addpath(genpath('functions'))
%% Define parameters
args.mu_0 = 4*pi*1e-7;
args.M    = 1.2706/args.mu_0;              % Magnetization   [A/m]
% permanent magnets
args.pm.L = 0.004;
args.pm.D = 0.002;
args.pm.m = args.M * (pi*args.pm.D^2/4*args.pm.L);
% micro robots (agents)
args.mr.L = 0.0004;
args.mr.D = 0.0002;
args.mr.m = args.M * (pi*args.mr.D^2/4*args.mr.L);
args.mr.m = args.mr.m * 1e3; %%%%%%%%% note
args.pm.L = args.pm.L * 1e0; %%%%%%%%% note
args.pm.D = args.pm.D * 1e0; %%%%%%%%% note
%
%
% values --> [x y z m_norm mu_0 m_agents]
a = [0.25 0.25 0.25];
z = 0.0;
phi = [30 150 270 -30 90 210]*pi/180;%phi = [30 150 270]*pi/180;
% phi = [000 180]*pi/180;
% phi = [000 090 180 270]*pi/180;
values = [
    a(1)*cos(phi(1)) a(1)*sin(phi(1)) z args.pm.m args.mu_0 args.mr.m
    a(2)*cos(phi(2)) a(2)*sin(phi(2)) z args.pm.m args.mu_0 args.mr.m
    a(3)*cos(phi(3)) a(3)*sin(phi(3)) z args.pm.m args.mu_0 args.mr.m
    a(1)*cos(phi(4)) a(1)*sin(phi(4)) z args.pm.m args.mu_0 args.mr.m
    a(2)*cos(phi(5)) a(2)*sin(phi(5)) z args.pm.m args.mu_0 args.mr.m
    a(3)*cos(phi(6)) a(3)*sin(phi(6)) z args.pm.m args.mu_0 args.mr.m
];
MagPos = [values(:,1:3) ones(size(values,1),1) zeros(size(values,1),2)];
for i=1:size(MagPos, 1)
    MagPos(i,4:6) = round(MagPos(i,4:6) ./ norm(MagPos(i,4:6)), 5);
end
vars = args;
vars.args = args;
vars.values = values;
vars.MagPos = MagPos;
vars.phi = phi;



%% Define parameters
domain = max(a) - 0.12;
step   = 0.01;
plotDomain = max(a) + 0.02;
spaceRegion = -domain:step:domain;
x_space = spaceRegion; z_space = x_space; y_space = x_space;
Npoints_space = length(x_space);
vars.x_space = x_space;
vars.y_space = y_space;
vars.z_space = z_space;
vars.plotDomain = plotDomain;
%
tspan = 0:1:500;%tspan = 0:1/5:120;%tspan = 0:1/24:30;
x_mr_ = min(x_space)+0.02:0.1:max(x_space)-0.02; % x_mr_ = [ 0.06 0.05 0.04 0.061 0.051 0.041];
y_mr_ = min(y_space)+0.02:0.1:max(y_space)-0.02; % y_mr_ = [ 0.06 0.05 0.04 0.061 0.051 0.041];
% x_mr_ = [ -0.010 -0.005 0.00 0.005 0.010 ];
% y_mr_ = [ -0.010 -0.005 0.00 0.005 0.010 ];
[x_mr__,y_mr__] = meshgrid(x_mr_,y_mr_);
x_mr_0 = reshape(x_mr__, 1, []);
y_mr_0 = reshape(y_mr__, 1, []);
%
vars.tspan = tspan;
vars.x_mr_0 = x_mr_0;
vars.y_mr_0 = y_mr_0;
%
counterName = 1;
cd('data')
fileName1 = '';
while strcmp(fileName1, '')
    if isfile([num2str(counterName) '_datas.mat'])
        counterName = counterName + 1;
    else
        fileName1 = [num2str(counterName) '_input.mat'];
        fileName2 = [num2str(counterName) '_movie.gif'];
        fileName3 = [num2str(counterName) '_datas.mat'];
    end
end
vars.fileName1 = fileName1;
vars.fileName2 = fileName2;
vars.fileName3 = fileName3;
cd('..')



dlgTitle    = 'PreCalculations';
dlgQuestion = 'Do you want symbolic calculations to be done ?';
choice = questdlg(dlgQuestion,dlgTitle,'Yes','No', 'Yes');
if strcmp(choice, 'Yes') == 1
    DoPreCalc = 1;
    symbolic_fun_b(values, fileName3, vars);
else
    DoPreCalc = 0;
end
vars.DoPreCalc = DoPreCalc;
% makeFunctionsFromSymbolic(b, b_m, str_psai, str_psai_sym) % inside _datas


%% calculate Psai_0 in order to have equilibrium in speceific point (using three or two magnets)
eq_point1 = [-0.05 0.05];
eq_point2 = [+0.03 +0.02];

% eq_point1 = [-0.06 0.05];
% Psai1 = 40*pi/180;
% Point = [-0.00 0.01];
% %%%
% % Psai_0 = calcPsai2(Point(1), Point(2), MagPos);
% %%%
% Psai_0 = calcPsai3(Point(1), Point(2), MagPos, Psai1);
% %%%
% 
% %%%%%%%%%%%%%%%
% [Psai_0,values] = calcPsaiNew(Point(1), Point(2), values, a, phi, z);
% MagPos = [values(:,1:3) zeros(size(values,1),3)];
% for i=1:size(MagPos, 1)
%     MagPos(i,4:6) = round(MagPos(i,4:6) ./ norm(MagPos(i,4:6)), 5);
% end
% %%%%%%%%%%%%%%%

% [Psai_0] = calculatePsai_2PM(eq_point1, MagPos);
% [Psai_0] = calculatePsai_3PM(eq_point1, MagPos);
% [Psai_0] = calculatePsai_4PM(eq_point1, eq_point2, MagPos);
[Psai_0] = calculatePsai_6PM(eq_point1, eq_point2, MagPos);


% calculateHesian_3PM(eq_point1, MagPos)

% %%%%%%%%%%%%%%%
% Psai_0 = zeros(size(values,1),1);
% Psai_0 = [
% %     30
%     150
%     -90
% ]*pi/180;
% Psai_0 = Psai_0 + [
% %     0
%     0
%     -75
% ]*pi/180;

% Psai_0 = 20*(pi/180) * ones(size(values,1),1);
% Psai_0 = (pi/180) * [20; 20; 60];
vars.Psai_0 = Psai_0;



%%
% printFig(vars, Psai_0, 0, [num2str(counterName) 'a']);
% printFig(vars, [0,0,0]*pi/180, 0);
printFig(vars, Psai_0, 0);



[Fx1, Fy1, Fz1] = CylFfield3(Psai_0,eq_point1);
[Bx1, By1, Bz1] = CylBfield3(Psai_0,eq_point1);

[Fx2, Fy2, Fz2] = CylFfield3(Psai_0,eq_point2);
[Bx2, By2, Bz2] = CylBfield3(Psai_0,eq_point2);


% [eq_x,eq_y] = findEqPoints_Dynamics(x_space,y_space,Psai_0);
[eq_x,eq_y] = findEqPoints_Minimization(x_space,y_space,Psai_0);
vars.eqPoints.x_eqPoints_0 = eq_x;
vars.eqPoints.y_eqPoints_0 = eq_y;



vars.dPath.indexInsideEqPoints = 1;
vars.dPath.dt = 0.2;
vars.dPath.startTime = 30;
vars.dPath.endTime = 60;
vars.dPath.maxDisp = 0.1;



idx = vars.dPath.indexInsideEqPoints;
% point_0 = [eq_x(idx); eq_y(idx)];
point_0 = [0; 0];
[x_desiredPath1,y_desiredPath1] = defineDesiredPath(point_0,vars.dPath);
vars.dPath.x_desiredPath1 = x_desiredPath1;
vars.dPath.y_desiredPath1 = y_desiredPath1;
point_0 = [-0.05; -0.05];
[x_desiredPath2,y_desiredPath2] = defineDesiredPath(point_0,vars.dPath);
vars.dPath.x_desiredPath2 = x_desiredPath2;
vars.dPath.y_desiredPath2 = y_desiredPath2;

[x_eqPoints_0] = designPsaiController(vars);

%%

cd('data')
save(fileName1, '-struct', 'vars')
cd('..')



% run_4__4
run_4__3___solveDynamicsAfterDefiningPathAndPsaiCont
if ~DoPreCalc
%     run_4__3
end




%%





function [Bx, By, Bz] = CylBfield3(Psai,Point)
    B_ = magnetic_field_symbolic(Point(1), Point(2), Psai);
    Bx = B_(1);
    By = B_(2);
    Bz = B_(3);
end
%
function [Fx, Fy, Fz] = CylFfield3(Psai,Point)
    F = force_field_symbolic(Point(1), Point(2), Psai);
    Fx = F(1);
    Fy = F(2);
    Fz = F(3);
end
%
function [MagPos] = sphericalToCartesian(MagPosR)

    for i_=1:size(MagPosR,1)
        r = MagPosR(i_,1);
        theta = MagPosR(i_,2);
        phi = MagPosR(i_,3);
        psai = MagPosR(i_,4);

        x_ = r*sin(theta)*cos(phi);
        y_ = r*sin(theta)*sin(phi);
        z_ = r*cos(theta);
        u_ = [x_ y_ z_]/r;

        x__ = -cos(theta)*cos(phi);
        y__ = -cos(theta)*sin(phi);
        z__ = sin(theta);
        u__ = [x__ y__ z__];

        R_mat = [
            cos(psai)+u_(1)^2*(1-cos(psai))  u_(1)*u_(2)*(1-cos(psai))-u_(3)*sin(psai) u_(1)*u_(3)*(1-cos(psai))+u_(2)*sin(psai);
            u_(1)*u_(2)*(1-cos(psai))+u_(3)*sin(psai)  cos(psai)+u_(2)^2*(1-cos(psai)) u_(2)*u_(3)*(1-cos(psai))-u_(1)*sin(psai);
            u_(1)*u_(3)*(1-cos(psai))-u_(2)*sin(psai) u_(2)*u_(3)*(1-cos(psai))+u_(1)*sin(psai) cos(psai)+u_(3)^2*(1-cos(psai));
        ];

        MagPos(i_,:) = [x_ y_ z_ (R_mat*u__')'];
    end

end
%
function [Psai] = calcPsai2(x, y, MagPos)

    for i=1:size(MagPos,1)
        x_ = MagPos(i,1);
        y_ = MagPos(i,2);
        z_ = MagPos(i,3);
        u_ = [x_ y_ z_];

        r_(i,:) = [x-x_ y-y_];
        r__(i,1) = sqrt( r_(i,1)^2 + r_(i,2)^2 + (0-z_)^2 );
    end
    c = -r__(1)^5/r__(2)^5;
    myFun = @(psai2)asin(c*sin(psai2))-acos(c*cos(psai2));
    psai2_0 = 40*pi/180;
    [psai2,fval,exitflag,output] = fsolve(myFun,psai2_0);
    [psai2,fval,exitflag,output] = fminsearch(@(psai2) abs(myFun(psai2)),psai2_0);
    psai1 = asin(c*sin(psai2));
    psai1 = acos(c*cos(psai2));
    error = asin(c*sin(psai2))-acos(c*cos(psai2));
    Psai= [psai1; psai2];
    %
    %
    c1 = 1/r__(1)^5;
    c2 = 1/r__(2)^5;
    myFun2 = @(psai)[c1*sin(psai(1))+c2*sin(psai(2)) c1*cos(psai(1))+c2*cos(psai(2))];
    psai_0 = [30 01]*pi/180;
    [psai12,fval,exitflag,output] = fsolve(myFun2,psai_0);
%     [psai2,fval,exitflag,output] = fminsearch(@(psai) abs(myFun2(psai)),psai_0);
    error = myFun2(psai12);
    Psai = psai12;
    %
    %
end
%
function [Psai] = calcPsai3(x, y, MagPos, Psai1)

    for i=1:size(MagPos,1)
        x_ = MagPos(i,1);
        y_ = MagPos(i,2);
        z_ = MagPos(i,3);
        u_ = [x_ y_ z_];

        r_(i,:) = [x-x_ y-y_];
        r__(i,1) = sqrt( r_(i,1)^2 + r_(i,2)^2 + (0-z_)^2 );
    end
    c32 = -r__(3)^5/r__(2)^5;
    c31 = -r__(3)^5/r__(1)^5;
    myFun = @(psai2)asin(c31*sin(Psai1)+c32*sin(psai2))-acos(c31*cos(Psai1)+c32*cos(psai2));
    psai2_0 = 30*pi/180;
    [psai2,fval,exitflag,output] = fsolve(myFun,psai2_0);
%     [psai2,fval,exitflag,output] = fminsearch(@(psai2) abs(myFun(psai2)),psai2_0);
    psai3 = asin(c31*sin(Psai1)+c32*sin(psai2));
    error = myFun(psai2);
    Psai = [Psai1; psai2; psai3];
    %
    %
    myFun2 = @(psai)asin(c31*sin(psai(1))+c32*sin(psai(2)))-acos(c31*cos(psai(1))+c32*cos(psai(2)));
    psai2_0 = [00; 00]*pi/180;
    [psai12,fval,exitflag,output] = fsolve(myFun2,psai2_0);
%     [psai12,fval,exitflag,output] = fminsearch(@(psai) abs(myFun2(psai)),psai2_0);
    psai3 = asin(c31*sin(psai12(1))+c32*sin(psai12(2)));
    Psai_ = [psai12(1); psai12(2); psai3];
    %
    %
    c1 = 1/r__(1)^5;
    c2 = 1/r__(2)^5;
    c3 = 1/r__(3)^5;
    myFun3 = @(psai)[c1*sin(psai(1))+c2*sin(psai(2))+c3*sin(psai(3)) c1*cos(psai(1))+c2*cos(psai(2))+c3*cos(psai(3))];
    psai_0 = [30 01 181]*pi/180;
    [psai123,fval,exitflag,output] = fsolve(myFun3,psai_0);
%     [psai2,fval,exitflag,output] = fminsearch(@(psai) abs(myFun3(psai)),psai_0);
    error = myFun3(psai123);
    Psai = psai123;
    %
    %
    %

    
end
%
function [Psai_0,values] = calcPsaiNew(x, y, values, a, phi, z)
    input_0 = [a 0.1 0.1 0.1];
    input_0 = [a(1) 0.1 0.1 0.1];
    opts = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt');
    [ans,fval,exitflag,output] = fsolve(@(input)myFun0(input, phi, x, y, z),input_0,opts);
    a = ans(1:3)';
    Psai_0 = ans(4:6)';
    values_ = [
        a(1)*cos(phi(1)) a(1)*sin(phi(1))
        a(2)*cos(phi(2)) a(2)*sin(phi(2))
        a(3)*cos(phi(3)) a(3)*sin(phi(3))
    ];
    values = [values_ values(:,3:end)];
end
%
function [ans] = myFun0(input, phi, x, y, z)
%     a(1) = input(1);
%     a(2) = input(2);
%     a(3) = input(3);
    a = [input(1) 0.25 0.25];
    psai(1) = input(2);
    psai(2) = input(3);
    psai(3) = input(4);
    for i=1:3
        x_ = a(i)*cos(phi(i));
        y_ = a(i)*sin(phi(i));
        z_ = z;

        r_(i,:) = [x-x_ y-y_];
        r__(i,1) = sqrt( r_(i,1)^2 + r_(i,2)^2 + (0-z_)^2 );
    end
    c1 = 1/r__(1)^5;
    c2 = 1/r__(2)^5;
    c3 = 1/r__(3)^5;
    ans = [
        c1*sin(psai(1))+c2*sin(psai(2))+c3*sin(psai(3))
        c1*cos(psai(1))+c2*cos(psai(2))+c3*cos(psai(3))
        ];
end
