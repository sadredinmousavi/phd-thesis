%% Defining parameters
 % all paramers have to be set in here
 % s 

clc, close all, clear all
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
cd('..');
addpath(genpath('functions'));
addpath(genpath('inputFiles'));
addpath(genpath('data'));


%% Define parameters
args.mu_0 = 4*pi*1e-7;
args.M    = 868 * 1e3;              % Magnetization   [A/m]   hcb n35
% permanent magnets
args.pm.L = 20 * 1e-3; %[m]
args.pm.D = 20 * 1e-3; %[m]
args.pm.V = args.pm.L ^ 3;
args.pm.m = args.M * args.pm.V;
% micro robots (agents)
args.mr.L = 0.5 * 1e-3; %[m]
args.mr.D = 0.5 * 1e-3; %[m]
args.mr.V = args.mr.L ^ 3; %(pi*args.mr.D^2/4*args.mr.L);
args.mr.m = args.M * args.mr.V;

args.mr.rho  = 7.5 * 1e3; %Rho_pla = 1.25 g/cm3 Rho_Neodymium = 7.5 g/cm3
args.mr.mass = args.mr.V * (args.mr.rho);
% free particles
args.fp.D = 0.010;
args.fp.rho  = 1.25 * 1e3; %Rho_pla = 1.25 g/cm3 Rho_Neodymium = 7.5 g/cm3
args.fp.mass = 4/3*pi*(args.fp.D/2)^3*(args.fp.rho);
% Lennard Jones potential
sigma = 1 * 1e-5;
% epsilun = 1e-10; % 0.01; 1e-8;
% threshold = 5*sigma;
epsilun = 10; % 1000; 1e-8;
threshold = 5*sigma + 0.002; %5*sigma + 0.004;
% epsilun = 5e-10; % 1000; 1e-8;
% threshold = 5*sigma; %5*sigma + 0.004;
% Drag
drag_coeff = 16/3;
viscosity = 0.01 * 5;


args.mr.m = args.mr.m * 1e0; %%%%%%%%% note
args.pm.L = args.pm.L * 1e0; %%%%%%%%% note
args.pm.D = args.pm.D * 1e0; %%%%%%%%% note
args.mr.mass = args.mr.mass * 1e0; %%%%%%%%% note
args.fp.mass = args.fp.mass * 1e0; %%%%%%%%% note



%% Define PMs locations
% values --> [x y z m_norm mu_0 m_agents]  for PMs
magNum = 4;
a = 0.14;
z = 0.0;
for i=1:magNum
    phi(i) = (0 + (i-1)*(360/magNum) )*(pi/180);
    values(i,:) = [a*cos(phi(i)) a*sin(phi(i)) z args.pm.m args.mu_0 args.mr.m];
end
% % calcPsaiFromEqFunc = @check_6PM; %check_6PM %check_6PM %check_6PM3 %calculatePsai_6PM
calcPsaiFromEqFunc = @calculatePsai_minimization; %@calculatePsai_4PM; %@calculatePsai_minimization;
symbolicFunctionHandle = @symbolic_fun_b;
findEqFromMinimization = 0;


%% Define space and plot parameters and options
domain = 0.15; % domain = max(a) - 0.12;
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
plotOptionsDyn.particlesTrack = 0;
plotOptionsDyn.channelLines = 1;
plotOptionsDyn.eqPoints = 0;
plotOptionsDyn.eqPointsTrack = 0;
plotOptionsDyn.printPsaiValues = 1;
plotOptionsDyn.printLambdaValues = 1;



%% Define MRs locations and parameters
r_mr = args.mr.D/2;
m_mr = args.mr.mass;
%
center1 = [-0.03, +0.03];
num1 = 3; %4
offset1 = 0.02;
x_serie1 = linspace(center1(1)-offset1, center1(1)+offset1, num1);
y_serie1 = linspace(center1(2)-offset1, center1(2)+offset1, num1);
[x_mr__,y_mr__] = meshgrid(x_serie1,y_serie1);
x_mr_1 = reshape(x_mr__, 1, []);
y_mr_1 = reshape(y_mr__, 1, []);
center2 = [+0.03, +0.03];
num2 = 3; %4
offset2 = 0.02;
x_serie2 = linspace(center2(1)-offset2, center2(1)+offset2, num2);
y_serie2 = linspace(center2(2)-offset2, center2(2)+offset2, num2);
[x_mr__,y_mr__] = meshgrid(x_serie2,y_serie2);
x_mr_2 = reshape(x_mr__, 1, []);
y_mr_2 = reshape(y_mr__, 1, []);

x_mr_0 = [x_mr_1 x_mr_2];
y_mr_0 = [y_mr_1 y_mr_2];
r_mr_0 = r_mr * ones(1, length(x_mr_0));
m_mr_0 = m_mr * ones(1, length(x_mr_0));
k_mr_0 = (16/3 * r_mr) * ones(2, length(x_mr_0));

%% Define Free particles locations
% % type 1 --> circle
% % type 2 --> rectngle
% % type 3 --> linear
r_fp = args.fp.D/2;
m_fp = args.fp.mass;
%


% x_fp_0 = [+0.00 +0.00];
% y_fp_0 = [+0.08 +0.04];
% t_fp_0 = [+0.00 +0.00];
% m_fp_0 = m_mr * ones(1, length(x_fp_0));
% i_fp_0 = m_mr * ones(1, length(x_fp_0));
% 
% fps{1}.type = 1;
% fps{1}.radius = 0.002;
% fps{2}.type = 2;
% fps{2}.points = 0.004 * [ [-1;-1] [+1;-1] [+1;+1] [-1;+1] ];


%

% x_fp_0 = [+0.00];
% y_fp_0 = [+0.05 ];
% t_fp_0 = [+0.00];
% m_fp_0 = m_fp * ones(1, length(x_fp_0));
% i_fp_0 = m_fp * ones(1, length(x_fp_0));
% %
% fps{1}.type = 1;
% fps{1}.radius = r_fp ;
% fps{1}.K = [16/3 16/3] * r_fp ;


% x_fp_0 = [+0.00];
% y_fp_0 = [+0.05];
% t_fp_0 = [+0.00];
% m_fp_0 = m_mr * ones(1, length(x_fp_0));
% i_fp_0 = m_mr * ones(1, length(x_fp_0));
% 
% fps{1}.type = 2;
% fps{1}.length = 0.010;
% fps{1}.K = [16/3 16/3] ;



% x_fp_0 = [+0.00];
% y_fp_0 = [+0.05];
% t_fp_0 = [+0.00];
% m_fp_0 = m_mr * ones(1, length(x_fp_0));
% i_fp_0 = m_mr * ones(1, length(x_fp_0));
% 
% fps{1}.type = 3;
% fps{1}.points = 0.004 * [ [-1;-1] [+1;-1] [+1;+1] [-1;+1] ];
% fps{1}.K = [16/3 16/3] ;



%% Define path and eqPoints desired locations during simulation time

% path inputs --> [point1 point2 startTime endTime dt]
% % eqPoint = [t1 t2 t3 ...; x1 x2 x3 ...; y1 y2 y3 ...]

% eqPoint1 = [ [0;-0.05;+0.00] [100;-0.05;+0.02] [200;-0.05;+0.04] [250;-0.03;+0.07] [300;-0.01;+0.09] [400;+0.00;+0.10] [500;+0.00;+0.05] [600;+0.00;+0.00] ];
% eqPoint2 = [ [0;+0.05;+0.00] [100;+0.05;+0.02] [200;+0.05;+0.04] [250;+0.03;+0.07] [300;+0.01;+0.09] [400;+0.00;+0.10] [500;+0.00;+0.05] [600;+0.00;+0.00] ];

% eqPoint1 = [ [0;-0.05;+0.08] [50;-0.05;+0.07] [100;-0.03;+0.06] [150;-0.00;+0.05] [200;-0.00;+0.03] [300;+0.00;+0.01] [400;+0.00;-0.01] [500;+0.00;-0.05] [600;+0.00;-0.07] [700;+0.03;-0.09] [800;+0.05;-0.09] [900;+0.09;-0.09] [1000;+0.13;-0.09] [1100;+0.14;-0.06] [1200;+0.14;-0.00] ];
% eqPoint2 = [ [0;+0.05;+0.08] [50;+0.05;+0.07] [100;+0.03;+0.06] [150;+0.00;+0.05] [200;+0.00;+0.03] [300;+0.00;+0.01] [400;+0.00;-0.01] [500;+0.00;-0.05] [600;+0.00;-0.07] [700;+0.03;-0.09] [800;+0.05;-0.09] [900;+0.09;-0.09] [1000;+0.13;-0.09] [1100;+0.14;-0.06] [1200;+0.14;-0.00] ];

% eqPoint1 = [ [0;-0.05;+0.08] [50;-0.05;+0.07] [70;-0.03;+0.06] [100;-0.00;+0.05] [130;-0.00;+0.02] [160;+0.00;-0.04] [200;+0.00;-0.07] [250;+0.03;-0.09] [300;+0.05;-0.09] [350;+0.09;-0.09] [400;+0.13;-0.09] [450;+0.14;-0.06] [500;+0.14;-0.00] ];
% eqPoint2 = [ [0;+0.05;+0.08] [50;+0.05;+0.07] [70;+0.03;+0.06] [100;+0.00;+0.05] [130;+0.00;+0.02] [160;+0.00;-0.04] [200;+0.00;-0.07] [250;+0.03;-0.09] [300;+0.05;-0.09] [350;+0.09;-0.09] [400;+0.13;-0.09] [450;+0.14;-0.06] [500;+0.14;-0.00] ];


% eqPoint1 = [ [0;-0.05;+0.00] [100;-0.05;+0.02] [200;-0.05;+0.04] [250;-0.03;+0.07] [300;-0.01;+0.09] [400;+0.00;+0.10] [500;+0.00;+0.08] ];
% eqPoint2 = [ [0;+0.05;+0.00] [100;+0.05;+0.02] [200;+0.05;+0.04] [250;+0.03;+0.07] [300;+0.01;+0.09] [400;+0.00;+0.10] [500;+0.00;+0.08] ];

eqPoint1 = [ [0;-0.00;+0.00] [100;-0.03;+0.00] [200;-0.05;+0.00] [200;-0.07;+0.00] [200;-0.05;+0.00] [100;-0.03;+0.00] [250;-0.00;-0.00] [250;-0.00;-0.03] [300;-0.00;-0.05] [400;-0.00;-0.07] [300;-0.00;-0.05] [250;-0.00;-0.03] [500;+0.00;+0.00] ];
eqPoint2 = [ [0;+0.00;+0.00] [100;+0.03;+0.00] [200;+0.05;+0.00] [200;+0.07;+0.00] [200;+0.05;+0.00] [100;+0.03;+0.00] [250;-0.00;-0.00] [250;+0.00;+0.03] [300;+0.00;+0.05] [400;+0.00;+0.07] [300;+0.00;+0.05] [250;+0.00;+0.03] [500;+0.00;+0.00] ];

% eqPoint1 = [ [0;-0.00;+0.00] [100;+0.02;+0.00] [200;+0.03;+0.02] [200;+0.03;+0.05] [250;-0.00;+0.05] [300;-0.03;+0.05] [400;-0.04;-0.00] [500;+0.00;+0.00] ];
% eqPoint2 = [ [0;+0.00;+0.00] [100;+0.02;+0.00] [200;+0.03;+0.02] [200;+0.03;+0.05] [250;+0.00;+0.05] [300;-0.03;+0.05] [400;-0.04;+0.00] [500;+0.00;+0.00] ];

% eqPoint1 = [ [0;-0.00;+0.00] [100;+0.02;+0.00] [200;+0.03;+0.02] [200;+0.03;+0.03] [250;-0.00;+0.03] [300;-0.03;+0.03] [400;-0.03;-0.00] [500;+0.00;+0.00] ];
% eqPoint2 = [ [0;+0.00;+0.00] [100;+0.02;+0.00] [200;+0.03;+0.02] [200;+0.03;+0.03] [250;+0.00;+0.03] [300;-0.03;+0.03] [400;-0.03;+0.00] [500;+0.00;+0.00] ];

eqPoints{1} = eqPoint1;
eqPoints{2} = eqPoint2;

%
%
endTime = max(eqPoint1(1,end),eqPoint2(1,end)) + 100;
startTime = 0;
stepForSolve = 0.01;
stepForOutput = 10;
tspan = startTime:stepForOutput:endTime;
%
%
usePreparedPsai = 0;%%%%% note


%% Define Walls

% path inputs --> [point1 point2 startTime endTime dt]
% wall = [type1 type2 type3 ...; x1 x2 x3 ...; y1 y2 y3 ...; x1_ x2_ x3_ ...; y1_ y2_ y3_ ...]
% wall_linear = [type; x1; y1; x2; y2]  --> type 0
% wall_circle = [type; x1; y1; r0; 00]  --> type 1

% walls_1 = [ [0;-0.04;+0.04;-0.03;-0.04] [0;-0.03;-0.04;-0.03;-0.09] [0;-0.03;-0.09;+0.09;-0.09] ];
% % walls_2 = [ [0;+0.04;+0.04;+0.03;-0.02] [0;+0.03;-0.02;+0.04;-0.03] [0;+0.04;-0.03;+0.09;-0.03] ];
% walls_2 = [ [0;+0.04;+0.04;+0.03;-0.04] [0;+0.03;-0.04;+0.09;-0.04] ];
% 
% walls = [ walls_1 walls_2];





%%
prepareTheInputs
solveTheSystem(vars);
plotOutputToGif(vars.fileName3);

% plotOptionsDyn.fieldVectors = 0;
% plotOptionsDyn.magnets = 1;
% plotOptionsDyn.areaBorders = 1;
% plotOptionsDyn.robotsTrack = 0;
% plotOptionsDyn.particlesTrack = 0;
% plotOptionsDyn.channelLines = 1;
% plotOptionsDyn.eqPoints = 0;
% plotOptionsDyn.eqPointsTrack = 0;
% plotOptionsDyn.printPsaiValues = 1;
% plotOptionsDyn.printLambdaValues = 1;
% plotOutputToGif(vars.fileName3, plotOptionsDyn);



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
