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
magNum = 10;
a = 0.25;
z = 0.0;
for i=1:magNum
    phi(i) = (0 + (i-1)*(360/magNum) )*(pi/180);
    values(i,:) = [a*cos(phi(i)) a*sin(phi(i)) z args.pm.m args.mu_0 args.mr.m];
end
% % calcPsaiFromEqFunc = @check_6PM; %check_6PM %check_6PM %check_6PM3 %calculatePsai_6PM
calcPsaiFromEqFunc = @calculatePsai_minimization;
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

% eqPoint1 = [ [000; +0.05; +0.00] [200; +0.05; +0.04] [400; +0.01; +0.09] [600; -0.03; +0.09] [800; -0.06; +0.04] [1000; -0.05; -0.00] ];% [1000; +0.05; +0.04]];
% eqPoint2 = [ [000; -0.05; -0.00] [200; -0.03; -0.06] [400; -0.02; -0.03] [600; +0.02; -0.03] [800; +0.02; -0.05] [1000; +0.05; -0.00] ];% [1000; -0.03; -0.06]];
% p1 = [1.39626340159546;1.39626340159546;0.779651245032826;1.57079632679490;1.12651856419339;1.57079632679490];
% p2 = [1.51940107369007;0;1.39124692092285;0.872664625997165;1.33007664670310;1.52982994555595];
% p3 = [1.22460040361176;0.872664625997165;1.58761411985301;0.349065850398866;1.54313506412893;1.25797162189381];
% p4 = [0.705451927114380;1.39626340159546;1.49705531561409;1.39626340159546;1.52434242494509;0.478744723660723];
% p5 = [6.28318530717959;1.58375123800553;1.16697743648885;1.50363231568715;5.41052068118242;0.356059539048953];
% Psai = [ p1 p2 p3 p4 p5 p1];

pointFunCyln2 = @(t,x,y,r,theta)[t; x+r*cos(theta*(pi/180)); y+r*sin(theta*(pi/180))];
pointNum = 16;
r = 0.03;
eqPoint1 = [ [000; +r; +0.00] ];
eqPoint2 = [ [000; -r; +0.00] ];
% theta = linspace(0,180,pointNum);
% for i=1:length(theta)
%     eqPoint1 = [eqPoint1 pointFunCyln2(eqPoint1(1,end)+50, 0,0,r, theta(i)) ];
%     eqPoint2 = [eqPoint2 pointFunCyln2(eqPoint2(1,end)+50, 0,0,r, theta(i)+180) ];
% end
eqPoint1 = [eqPoint1 pointFunCyln2(eqPoint1(1,end)+50, 0,0,r, 15) ];
eqPoint2 = [eqPoint2 pointFunCyln2(eqPoint2(1,end)+50, 0,0,r, 180-15) ];
%
eqPoint1 = [eqPoint1 pointFunCyln2(eqPoint1(1,end)+50, 0,0,r, 30) ];
eqPoint2 = [eqPoint2 pointFunCyln2(eqPoint2(1,end)+50, 0,0,r, 180-30) ];
%
eqPoint1 = [eqPoint1 pointFunCyln2(eqPoint1(1,end)+50, 0,0,r, 45) ];
eqPoint2 = [eqPoint2 pointFunCyln2(eqPoint2(1,end)+50, 0,0,r, 180-45) ];
%
eqPoint1 = [eqPoint1 pointFunCyln2(eqPoint1(1,end)+50, 0.02,0.02,r, 45) ];
eqPoint2 = [eqPoint2 pointFunCyln2(eqPoint2(1,end)+50, 0.02,0.02,r, 180-45) ];
%
eqPoint1 = [eqPoint1 pointFunCyln2(eqPoint1(1,end)+50, 0.04,0.04,r, 45) ];
eqPoint2 = [eqPoint2 pointFunCyln2(eqPoint2(1,end)+50, 0.04,0.04,r, 180-45) ];
%
% for i=1:length(theta)
%     eqPoint1 = [eqPoint1 pointFunCyln2(eqPoint1(1,end)+50, 0.05,0.05,r, theta(i)) ];
%     eqPoint2 = [eqPoint2 pointFunCyln2(eqPoint2(1,end)+50, 0.05,0.05,r, theta(i)+180) ];
% end
% eqPoint1 = [eqPoint1 pointFunCyln2(eqPoint1(1,end)+50, 0.02,0,r, theta(i)) ];
% eqPoint2 = [eqPoint2 pointFunCyln2(eqPoint2(1,end)+50, 0.02,0,r, theta(i)+180) ];
    

%
%
tspan = 0:10:max(eqPoint1(1,end),eqPoint2(1,end))+600;%tspan = 0:1/5:120;%tspan = 0:1/24:30;
%
%

usePrepaidPsai = 0;%%%%% note

eqPointsNum = 2;



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
