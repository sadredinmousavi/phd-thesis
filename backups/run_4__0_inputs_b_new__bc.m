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
vars = args;
vars.args = args;

%% Define PMs locations
% values --> [x y z m_norm mu_0 m_agents]  for PMs
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
%     a(1)*cos(phi(5)) a(1)*sin(phi(5)) z args.pm.m args.mu_0 args.mr.m
%     a(2)*cos(phi(6)) a(2)*sin(phi(6)) z args.pm.m args.mu_0 args.mr.m
];
MagPos = [values(:,1:3) ones(size(values,1),1) zeros(size(values,1),2)];
for i=1:size(MagPos, 1)
    MagPos(i,4:6) = round(MagPos(i,4:6) ./ norm(MagPos(i,4:6)), 5);
end
%
%
vars.values = values;
vars.MagPos = MagPos;
vars.phi = phi;
%
vars.calcPsaiFromEqFunc = @check_6PM; %check_6PM %check_6PM %check_6PM3 %calculatePsai_6PM
vars.findEqFromMinimization = 0;



%% Define filing configs
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
%
%
%
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




%% Define space and plot parameters and options
domain = max(a) - 0.12;
step   = 0.01;
plotDomain = max(a) + 0.02;
spaceRegion = -domain:step:domain;
x_space = spaceRegion; z_space = x_space; y_space = x_space;
Npoints_space = length(x_space);
%
%
%
vars.x_space = x_space;
vars.y_space = y_space;
vars.z_space = z_space;
vars.plotDomain = plotDomain;
%
vars.plotOptions.dynamic.fieldVectors = 0;
vars.plotOptions.dynamic.magnets = 1;
vars.plotOptions.dynamic.areaBorders = 1;
vars.plotOptions.dynamic.robotsTrack = 0;
vars.plotOptions.dynamic.eqPoints = 1;
vars.plotOptions.dynamic.eqPointsTrack = 1;



%% Define MRs locations
% x_mr_ = min(x_space)+0.02:0.1:max(x_space)-0.02; % x_mr_ = [ 0.06 0.05 0.04 0.061 0.051 0.041];
% y_mr_ = min(y_space)+0.02:0.1:max(y_space)-0.02; % y_mr_ = [ 0.06 0.05 0.04 0.061 0.051 0.041];
% [x_mr__,y_mr__] = meshgrid(x_mr_,y_mr_);
% x_mr_0 = reshape(x_mr__, 1, []);
% y_mr_0 = reshape(y_mr__, 1, []);
center1 = [-0.05, 0.02];
num1 = 2; %4
offset1 = 0.02;
x_serie1 = linspace(center1(1)-offset1, center1(1)+offset1, num1);
y_serie1 = linspace(center1(2)-offset1, center1(2)+offset1, num1);
[x_mr__,y_mr__] = meshgrid(x_serie1,y_serie1);
x_mr_1 = reshape(x_mr__, 1, []);
y_mr_1 = reshape(y_mr__, 1, []);
% center2 = [+0.05, 0.02];
% num2 = 2; %4
% offset2 = 0.02;
% x_serie2 = linspace(center2(1)-offset2, center2(1)+offset2, num2);
% y_serie2 = linspace(center2(2)-offset2, center2(2)+offset2, num2);
% [x_mr__,y_mr__] = meshgrid(x_serie2,y_serie2);
% x_mr_2 = reshape(x_mr__, 1, []);
% y_mr_2 = reshape(y_mr__, 1, []);
% x_mr_0 = [x_mr_1 x_mr_2];
% y_mr_0 = [y_mr_1 y_mr_2];
x_mr_0 = [x_mr_1 ];
y_mr_0 = [y_mr_1 ];
%
%
vars.x_mr_0 = x_mr_0;
vars.y_mr_0 = y_mr_0;


%% calculate Psai_0 in order to have equilibrium in 2 points (using 4 or 6 magnets)

dlgTitle    = 'Hand Calculations';
dlgQuestion = 'Do you want eqPoint calculations to be done ?';
choice = questdlg(dlgQuestion,dlgTitle,'Yes','No', 'Yes');
if strcmp(choice, 'Yes') == 1
    % eqPoint1 = [ [0;+0.05;+0.00] [300;+0.03;+0.04] [600;+0.03;+0.09] ];
    % eqPoint2 = [ [0;-0.11;-0.02] [300;-0.10;-0.09] [600;+0.07;-0.10] ];
    eq_point1 = [+0.03 +0.04];
    eq_point2 = [-0.10 -0.09];
    vars.plotOptions.static.plotEqPoints = 1;
    vars.plotOptions.static.eq_point1 = eq_point1;
    vars.plotOptions.static.eq_point2 = eq_point2;
    

%     [rankM, error, hasAns, isStable, Psai, hessian] = vars.calcPsaiFromEqFunc(eq_point1, eq_point2, MagPos);
    
    PsaiSerie =  [1.50083706544196,1.49616104582681,1.48226502196774,1.45953645336855,1.42859437596614,1.39025589550213,1.39025589550213,1.42859437596614,1.45953645336855,1.48226502196774,1.49616104582681,1.50083706544196;0,0,0,0,0,0,0,0,0,0,0,0;1.23454041983653,1.23420295400830,1.23320057514786,1.23156306013890,1.22933911331976,1.22659500225686,1.22659500225686,1.22933911331976,1.23156306013890,1.23320057514786,1.23420295400830,1.23454041983653;0,0.174532925199433,0.349065850398866,0.523598775598299,0.698131700797732,0.872664625997165,5.41052068118242,5.58505360638185,5.75958653158129,5.93411945678072,6.10865238198015,6.28318530717959;1.29613062459016,1.30396670402988,1.32714036785550,1.36470110043027,1.41522106858842,1.47696937203877,1.47696937203877,1.41522106858842,1.36470110043027,1.32714036785550,1.30396670402988,1.29613062459016;1.54587681274818,1.54617198511749,1.54704852101410,1.54847975264256,1.55042214401606,1.55281662500641,1.55281662500641,1.55042214401606,1.54847975264256,1.54704852101410,1.54617198511749,1.54587681274818];
    for cnt=1:length(PsaiSerie)
        Psai = PsaiSerie(:,cnt);
        printFig(vars, Psai, 0);
    end
    
    %
    [Fx1, Fy1, Fz1] = CylFfield3(Psai,eq_point1);
    [Bx1, By1, Bz1] = CylBfield3(Psai,eq_point1);
    %
    [Fx2, Fy2, Fz2] = CylFfield3(Psai,eq_point2);
    [Bx2, By2, Bz2] = CylBfield3(Psai,eq_point2);
    %
    printFig(vars, Psai, 0);
    %
    [V,D] = eig(hessian.point1); % V(:,i)
    angle1_point1 = atan(V(2,1)/V(1,1))*(180/pi);
    angle2_point1 = atan(V(2,2)/V(1,2))*(180/pi);
    d1_point1 = D(1,1);
    d2_point1 = D(2,2);
    [V,D] = eig(hessian.point2); % V(:,i)
    angle1_point2 = atan(V(2,1)/V(1,1))*(180/pi);
    angle2_point2 = atan(V(2,2)/V(1,2))*(180/pi);
    d1_point2 = D(1,1);
    d2_point2 = D(2,2);



    newStep = 0.01;
    figure
    hold on
    plot(MagPos(:,1), MagPos(:,2), 'r.', 'MarkerSize', 15)
    plot(eq_point1(1), eq_point1(2), '+r')
    for x=-domain:newStep:domain
        for y=-domain:newStep:domain
            [rankM, error, hasAns, isStable, Psai, hessian] = vars.calcPsaiFromEqFunc(eq_point1, [x y], MagPos);
    %         if hasAns
            if isStable && hasAns
                plot(x, y, 'co')
            else
                plot(x, y, 'ko')
            end
            if hasAns
                a=1;
            end
        end
    end
    axis square
end

%% calculate Psai_0 in order to have equilibrium in 3 points (using 6 magnets)
% eq_point1 = [-0.02 +0.04];%eq_point1 = [-0.01 +0.05];
% eq_point2 = [+0.04 -0.11];%eq_point2 = [+0.04 -0.13];
% eq_point3 = [+0.04 -0.11];%eq_point3 = [+0.06 -0.11];
% 
% [rankM, error, hasAns, isStable, Psai] = check_6PM3(eq_point1, eq_point2, eq_point3, MagPos);
% printFig(vars, Psai, 0);
% newStep = 0.01;
% 
% allPointCount = 0;
% equalPointCount = 0;
% for x2=-domain:newStep:domain
%     for y2=-domain:newStep:domain
%         x3=x2;y3=y2;
% %         for x3=-domain:newStep:domain
% %             for y3=-domain:newStep:domain
%                 [rankM,error, hasAns, isStable, Psai] = check_6PM3(eq_point1, [x2 y2], [x3 y3], MagPos);
% %                 if hasAns
%                 if isStable && hasAns
%                     allPointCount = allPointCount + 1;
%                     eqpoint1(allPointCount, :) = eq_point1;
%                     eqpoint2(allPointCount, :) = [x2 y2];
%                     eqpoint3(allPointCount, :) = [x3 y3];
%                     if x2==x3 && y2==y3
%                         equalPointCount = equalPointCount + 1;
%                     end
%                 end
%                 if hasAns
%                     a=1;
%                 end
% %             end
% %         end
%     end
% end
% %
% figure
% hold on
% plot(MagPos(:,1), MagPos(:,2), 'r.', 'MarkerSize', 15)
% plot(eq_point1(1), eq_point1(2), '+k')
% for x=-domain:newStep:domain
%     for y=-domain:newStep:domain
%         plot(x, y, 'ko')
%     end
% end
% for i=1:length(eqpoint2)
%     plot(eqpoint2(i,1), eqpoint2(i,2), 'co')
% end



%% Define path and eqPoints desired locations during simulation time

% path inputs --> [point1 point2 startTime endTime dt]

% % eqPoint = [t1 t2 t3 ...; x1 x2 x3 ...; y1 y2 y3 ...]
% eqPoint1 = [ [0;-0.05;0.05] [300;-0.02;0.02] ];
% eqPoint2 = [ [0;+0.05;0.05] [300;+0.02;0.02] ];

%
%

% p1 = [-0.05 0.05];
% p2 = [+0.05 0.05];
% [path1] = defineDesiredPath(p1, p2, 0, 10, 1);
% p1 = [+0.05 0.02];
% p2 = [-0.02 0.04];
% [path2] = defineDesiredPath(p1, p2, 0, 10, 1);
% eqPoint1 = [path1 [20;.1;.1] [30;.2;.2] ];
% eqPoint2 = [path2];

%
%

% r = 0.1;%0.5
% pointFun = @(t,theta)[t; r*cos(theta*(pi/180)); r*sin(theta*(pi/180))];
% % eqPoint1 = [ pointFun(0,0) pointFun(100,15) pointFun(200,30) pointFun(400,60)];
% % eqPoint2 = [ pointFun(0,180) pointFun(100,165) pointFun(200,150) pointFun(400,120)];
% deg = [0:15:75 60:-15:0]; %[0:15:90 75:-15:0];
% timeStep = 150;%100
% eqPoint1 = pointFun(0,deg(1));
% eqPoint2 = pointFun(0,180-deg(1));
% for i=1:length(deg)-1
%     eqPoint1 = [eqPoint1 pointFun(timeStep*i,deg(i+1))];
%     eqPoint2 = [eqPoint2 pointFun(timeStep*i,180-deg(i+1))];
% end
% eqPoint1 = [eqPoint1 [eqPoint1(1,end)+timeStep; -0.1; 0] [eqPoint1(1,end)+timeStep; -0.14; 0] ];
% eqPoint2 = [eqPoint2 [eqPoint2(1,end)+timeStep; +0.1; 0] [eqPoint2(1,end)+timeStep; +0.14; 0] ];

%
%

% pointFunCyln = @(t,r,theta)[t; r*cos(theta*(pi/180)); r*sin(theta*(pi/180))];
% pointFunCart = @(t,x,y)[t; x; y];
% speed = 0.10 * 1e-3; %m/s
% refreshTime = 50.0; %s
% %
% r = 0.1;
% deg0 = 0;
% deg1 = 70;
% t0 = 300;
% path1 = defineCirclePath(r, [deg0 deg1], t0, refreshTime, speed);
% path2 = defineCirclePath(r, 180-[deg0 deg1], t0, refreshTime, speed);
% eqPoint1 = [ [0; path1(2,1); path1(3,1)] path1];
% eqPoint2 = [ [0; path2(2,1); path2(3,1)] path2];
% 
% path1 = defineLinearPath(eqPoint1(2:3,end), eqPoint1(2:3,end)-[0;0.2], eqPoint1(1,end), refreshTime, speed);
% path2 = defineLinearPath(eqPoint2(2:3,end), eqPoint2(2:3,end)-[0;0.2], eqPoint2(1,end), refreshTime, speed);
% eqPoint1 = [ eqPoint1 path1];
% eqPoint2 = [ eqPoint2 path2];

%
%

% eqPoint1 = [eqPoint1 [eqPoint1(1,end)+timeStep; -0.1; 0] [eqPoint1(1,end)+timeStep; -0.14; 0] ];
% eqPoint2 = [eqPoint2 [eqPoint2(1,end)+timeStep; +0.1; 0] [eqPoint2(1,end)+timeStep; +0.14; 0] ];

%
%

% eqPoint1 = [ [0;+0.05;+0.00] [300;+0.03;+0.04] [600;+0.03;+0.09] ];
% eqPoint2 = [ [0;-0.11;-0.02] [300;-0.10;-0.09] [600;+0.07;-0.10] ];

%
%

% eqPoint1 = [ [0;-0.06;+0.00] [300;-0.06;+0.03] [600;+0.03;+0.09] ];
% eqPoint2 = [ [0;+0.10;-0.00] [300;+0.06;-0.03] [600;+0.07;-0.10] ];
% eqPoint3 = [ [0;+0.10;-0.00] [300;+0.06;-0.03] [600;+0.07;-0.10] ];

%
%

% eqPoint1 = [ [000; +0.00; 0] [200; +0.0; +0.05] [400; +0.1; +0.05] [600; +0.14; +0.05]];
% eqPoint2 = [ [000; -0.00; 0] [200; -0.0; +0.05] [400; -0.1; +0.05] [600; -0.14; +0.05]];

%
%

eqPoint1 = [ [000; +0.00; 0] [200; +0.0; +0.05] [400; +0.1; +0.05] [600; +0.14; +0.05]];
eqPoint2 = [ [000; -0.00; 0] [200; -0.0; +0.05] [400; -0.1; +0.05] [600; -0.14; +0.05]];

%
%
tspan = 0:10:max(eqPoint1(1,end),eqPoint2(1,end))+600;%tspan = 0:1/5:120;%tspan = 0:1/24:30;
%
%
vars.tspan = tspan;
vars.eqPointsNum = 2;
vars.eqPoint1 = eqPoint1;
vars.eqPoint2 = eqPoint2;


% % [eq_x,eq_y] = findEqPoints_Dynamics(x_space,y_space,Psai_0);
% [eq_x,eq_y] = findEqPoints_Minimization(x_space,y_space,Psai_0);
% vars.eqPoints.x_eqPoints_0 = eq_x;
% vars.eqPoints.y_eqPoints_0 = eq_y;

% 
% designPsaiController(vars);
% 


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
% [Psai_0] = calculatePsai_6PM(eq_point1, eq_point2, MagPos);


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
[rankM, error, hasAns, isStable, Psai_0] = vars.calcPsaiFromEqFunc(eq_point1, eq_point2, MagPos);
vars.Psai_0 = Psai_0;



%%
% printFig(vars, Psai_0, 0, [num2str(counterName) 'a']);
% printFig(vars, [0,0,0]*pi/180, 0);
% printFig(vars, Psai_0, 0);



[Fx1, Fy1, Fz1] = CylFfield3(Psai_0,eq_point1);
[Bx1, By1, Bz1] = CylBfield3(Psai_0,eq_point1);

[Fx2, Fy2, Fz2] = CylFfield3(Psai_0,eq_point2);
[Bx2, By2, Bz2] = CylBfield3(Psai_0,eq_point2);





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

function [eqPoint1] = defineLinearPath(eqPoint0, eqPoint1, t0, refreshTime, speed)
    point1_x0 = eqPoint0(1);
    point1_y0 = eqPoint0(2);
    point1_x1 = eqPoint1(1);
    point1_y1 = eqPoint1(2);
    %
    speed_x = speed * (point1_x1-point1_x0)/sqrt((point1_x1-point1_x0)^2+(point1_y1-point1_y0)^2);
    speed_y = speed * (point1_y1-point1_y0)/sqrt((point1_x1-point1_x0)^2+(point1_y1-point1_y0)^2);
    t1 = max( abs( (point1_x1-point1_x0) / speed_x), abs( (point1_y1-point1_y0) / speed_y) );
    time = 0:refreshTime:t1;
    eqPoint1 = zeros(3,length(time));
    for i=1:length(time)
        eqPoint1(:,i) = [t0+time(i); point1_x0+speed_x*time(i); point1_y0+speed_y*time(i)];
    end
end
%
function [eqPoint1] = defineCirclePath(r, deg, t0, refreshTime, speed)
    pointFunCyln = @(t,r,theta)[t; r*cos(theta*(pi/180)); r*sin(theta*(pi/180))];
    omega = speed/r*(180/pi); %deg/s
    t1 = abs( (deg(2)-deg(1)) / omega );
    time = 0:refreshTime:t1;
    eqPoint1 = zeros(3,length(time));
    for i=1:length(time)
        eqPoint1(:,i) = pointFunCyln(t0+time(i),r,deg(1)+omega*time(i));
    end
end
%
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
