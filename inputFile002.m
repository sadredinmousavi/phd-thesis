%% Defining parameters
 % all paramers have to be set in here
 % s 


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
%
%
calcPsaiFromEqFunc = @check_6PM; %check_6PM %check_6PM %check_6PM3 %calculatePsai_6PM
symbolicFunctionHandle = @symbolic_fun_b;
findEqFromMinimization = 0;



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
%
plotOptionsDyn.fieldVectors = 0;
plotOptionsDyn.magnets = 1;
plotOptionsDyn.areaBorders = 1;
plotOptionsDyn.robotsTrack = 0;
plotOptionsDyn.eqPoints = 1;
plotOptionsDyn.eqPointsTrack = 1;
plotOptionsDyn.printPsaiValues = 1;
plotOptionsDyn.printLambdaValues = 1;



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

usePrepaidPsai = 1;%%%%% note

eqPointsNum = 2;


