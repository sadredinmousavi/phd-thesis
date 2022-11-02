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
phi = [30 150 270 -30 90 210]*pi/180;
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

center1 = [-0.00, 0.00];
num1 = 4; %4
offset1 = 0.02;
x_serie1 = linspace(center1(1)-offset1, center1(1)+offset1, num1);
y_serie1 = linspace(center1(2)-offset1, center1(2)+offset1, num1);
[x_mr__,y_mr__] = meshgrid(x_serie1,y_serie1);
x_mr_0 = reshape(x_mr__, 1, []);
y_mr_0 = reshape(y_mr__, 1, []);





%% Define path and eqPoints desired locations during simulation time

% path inputs --> [point1 point2 startTime endTime dt]

% % eqPoint = [t1 t2 t3 ...; x1 x2 x3 ...; y1 y2 y3 ...]
% eqPoint1 = [ [0;-0.05;0.05] [300;-0.02;0.02] ];
% eqPoint2 = [ [0;+0.05;0.05] [300;+0.02;0.02] ];


% eqPoint1 = [ [000; +0.05; +0.00] [200; +0.05; +0.04] [400; +0.1; +0.05] [600; +0.14; +0.05]];
% eqPoint2 = [ [000; -0.05; +0.00] [200; -0.03; -0.06] [400; -0.1; +0.05] [600; -0.14; +0.05]];


% eqPoint1 = [ [000; +0.05; +0.04] [200; +0.05; +0.04] [400; +0.05; +0.04] [600; +0.05; +0.04] [800; +0.05; +0.04] [1000; +0.05; +0.04]];
% eqPoint2 = [ [000; -0.03; -0.06] [200; -0.03; -0.06] [400; -0.03; -0.06] [600; -0.03; -0.06] [800; -0.03; -0.06] [1000; -0.03; -0.06]];
% 
% Psai = [1.51940107369007,1.51940107369007,1.52199133023939,1.52199133023939,1.52968153257363,1.50322683299611;0,0,0.174532925199433,0.174532925199433,0.349065850398866,0.349065850398866;1.39124692092285,1.39124692092285,1.39500293380128,1.39500293380128,1.40614231089409,1.39154256778066;0.872664625997165,5.41052068118242,0.872664625997165,5.41052068118242,0.872664625997165,1.04719755119660;1.33007664670310,1.33007664670310,1.32691320467839,1.32691320467839,1.31750416290284,1.41746590014614;1.52982994555595,1.52982994555595,1.52717529388934,1.52717529388934,1.51929011301187,1.56719160308122];

eqPoint1 = [ [000; +0.05; +0.00] [200; +0.05; +0.04] [400; +0.01; +0.09] [600; -0.03; +0.09] [800; -0.06; +0.04] [1000; -0.05; -0.00] ];% [1000; +0.05; +0.04]];
eqPoint2 = [ [000; -0.05; -0.00] [200; -0.03; -0.06] [400; -0.02; -0.03] [600; +0.02; -0.03] [800; +0.02; -0.05] [1000; +0.05; -0.00] ];% [1000; -0.03; -0.06]];
p1 = [1.39626340159546;1.39626340159546;0.779651245032826;1.57079632679490;1.12651856419339;1.57079632679490];
p2 = [1.51940107369007;0;1.39124692092285;0.872664625997165;1.33007664670310;1.52982994555595];
p3 = [1.22460040361176;0.872664625997165;1.58761411985301;0.349065850398866;1.54313506412893;1.25797162189381];
p4 = [0.705451927114380;1.39626340159546;1.49705531561409;1.39626340159546;1.52434242494509;0.478744723660723];
p5 = [6.28318530717959;1.58375123800553;1.16697743648885;1.50363231568715;5.41052068118242;0.356059539048953];
Psai = [ p1 p2 p3 p4 p5 p1];

%
%
tspan = 0:10:max(eqPoint1(1,end),eqPoint2(1,end))+600;%tspan = 0:1/5:120;%tspan = 0:1/24:30;
%
%

usePrepaidPsai = 1;%%%%% note

eqPointsNum = 2;


