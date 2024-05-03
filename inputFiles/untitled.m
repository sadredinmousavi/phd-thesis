%% Defining parameters
 % all paramers have to be set in here
 % s 

clc, close all, clear all
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
cd('..');
addpath(genpath('class'));


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
args.mr.k = (16/3 * args.mr.D/2);
args.mr.rho  = 7.5 * 1e3; %Rho_pla = 1.25 g/cm3 Rho_Neodymium = 7.5 g/cm3
args.mr.mass = args.mr.V * (args.mr.rho);

% free particles
args.fp.rho  = 1.25 * 1e3; %Rho_pla = 1.25 g/cm3 Rho_Neodymium = 7.5 g/cm3

% Lennard Jones potential
args.sigma = 1 * 1e-5;
% args.epsilun = 1e-10; % 0.01; 1e-8;
% args.threshold = 5*sigma;
args.epsilun = 10; % 1000; 1e-8;
args.threshold = 5*args.sigma + 0.002; %5*sigma + 0.004;
% args.epsilun = 5e-10; % 1000; 1e-8;
% args.threshold = 5*sigma; %5*sigma + 0.004;

% Drag
args.drag_coeff = 16/3;
args.viscosity = 0.1 * 1;

% manipulation
args.mr.m = args.mr.m * 1e0; %%%%%%%%% note
args.pm.L = args.pm.L * 1e0; %%%%%%%%% note
args.pm.D = args.pm.D * 1e0; %%%%%%%%% note
args.mr.mass = args.mr.mass * 1e0; %%%%%%%%% note
args.fp.mass = args.fp.mass * 1e0; %%%%%%%%% note



%% Define PMs locations
% values --> [x y z m_norm mu_0 m_agents]  for PMs
magNum = 8;
a = 0.25;%14
z = 0.0;
for i=1:magNum
    phi(i) = (0 + (i-1)*(360/magNum) )*(pi/180);
    psai(i) = 0;
    MagPos(i,:) = [a*cos(phi(i)) a*sin(phi(i)) z 1 0 0];
%     MagPos(i,4:6) = round(MagPos(i,4:6) ./ norm(MagPos(i,4:6)), 5);
end

configs.epms.magNum = magNum;
configs.epms.a = a;
configs.epms.z = z;
configs.epms.MagPos = MagPos;
configs.epms.phi = phi;
configs.epms.psai = psai;



%% Define space and plot parameters and options
domain = 0.25; % domain = max(a) - 0.12;
step   = 0.01;
plotDomain = max(a) + 0.02;
spaceRegion = -domain:step:domain;
x_space = spaceRegion; z_space = x_space; y_space = x_space;
Npoints_space = length(x_space);


configs.vars.domain = domain;
configs.vars.step = step;
configs.vars.plotDomain = plotDomain;
configs.vars.Npoints_space = Npoints_space;
configs.vars.x_space = x_space;
configs.vars.y_space = y_space;
configs.vars.z_space = z_space;


%% Plot Options
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


plotOptions.dynamic = plotOptionsDyn;



%% Define functions
% % calcPsaiFromEqFunc = @check_6PM; %check_6PM %check_6PM %check_6PM3 %calculatePsai_6PM
handles.calcPsaiFromEqFunc = @calculatePsai_minimization;%@calculatePsai_6PM; %@calculatePsai_4PM; %@calculatePsai_minimization;
handles.symbolicFunctionHandle = @symbolic_fun_b;
handles.findEqFromMinimization = 1;


%% Define MRs locations and parameters
groups{1}.center = [-0.03, +0.03];
groups{1}.num = 3; %4
groups{1}.offset = 0.02;
%
groups{2}.center = [+0.03, +0.03];
groups{2}.num = 3; %4
groups{2}.offset = 0.02;

configs.mrs.seedData = groups;


%% Define path and eqPoints desired locations during simulation time

% path inputs --> [point1 point2 startTime endTime dt]
% % eqPoint = [t1 t2 t3 ...; x1 x2 x3 ...; y1 y2 y3 ...]

eqPoint1 = [ [0;-0.01;+0.01] [10;-0.03;+0.00] [20;-0.05;+0.00] [30;-0.07;+0.00] [40;-0.05;+0.00] [50;-0.03;+0.00] [60;-0.00;-0.00] [70;-0.00;-0.03] [80;-0.00;-0.05] [90;-0.00;-0.07] [100;-0.00;-0.05] [120;-0.00;-0.03] [130;+0.00;+0.00] ];
eqPoint2 = [ [0;+0.01;+0.01] [10;+0.03;+0.00] [20;+0.05;+0.00] [30;+0.07;+0.00] [40;+0.05;+0.00] [50;+0.03;+0.00] [60;-0.00;-0.00] [70;+0.00;+0.03] [80;+0.00;+0.05] [90;+0.00;+0.07] [100;+0.00;+0.05] [120;+0.00;+0.03] [130;+0.00;+0.00] ];

% eqPoint1 = [ [0;-0.00;+0.00] [100;+0.02;+0.00] [200;+0.03;+0.02] [200;+0.03;+0.05] [250;-0.00;+0.05] [300;-0.03;+0.05] [400;-0.04;-0.00] [500;+0.00;+0.00] ];
% eqPoint2 = [ [0;+0.00;+0.00] [100;+0.02;+0.00] [200;+0.03;+0.02] [200;+0.03;+0.05] [250;+0.00;+0.05] [300;-0.03;+0.05] [400;-0.04;+0.00] [500;+0.00;+0.00] ];

rr = 0.03;
eqPoint1 = [ [0;-0.00;+0.00] [10;rr*cosd(000+000);rr*sind(000+000)] [20;rr*cosd(030+000);rr*sind(030+000)] [30;rr*cosd(060+000);rr*sind(060+000)] [40;rr*cosd(90+000);rr*sind(090+000)] [50;rr*cosd(120+000);rr*sind(120+000)] [60;rr*cosd(150+000);rr*sind(150+000)] [70;rr*cosd(180+000);rr*sind(180+000)] [80;rr*cosd(210+000);rr*sind(210+000)] [90;rr*cosd(240+000);rr*sind(240+000)] [100;rr*cosd(270+000);rr*sind(270+000)] [110;rr*cosd(300+000);rr*sind(300+000)] [120;rr*cosd(330+000);rr*sind(330+000)] [130;rr*cosd(360+000);rr*sind(360+000)] [140;-0.00;+0.00] ];
eqPoint2 = [ [0;-0.00;+0.00] [10;rr*cosd(000+000);rr*sind(000+000)] [20;rr*cosd(030+000);rr*sind(030+000)] [30;rr*cosd(060+000);rr*sind(060+000)] [40;rr*cosd(90+000);rr*sind(090+000)] [50;rr*cosd(120+000);rr*sind(120+000)] [60;rr*cosd(150+000);rr*sind(150+000)] [70;rr*cosd(180+000);rr*sind(180+000)] [80;rr*cosd(210+000);rr*sind(210+000)] [90;rr*cosd(240+000);rr*sind(240+000)] [100;rr*cosd(270+000);rr*sind(270+000)] [110;rr*cosd(300+000);rr*sind(300+000)] [120;rr*cosd(330+000);rr*sind(330+000)] [130;rr*cosd(360+000);rr*sind(360+000)] [140;-0.00;+0.00] ];

eqPoints{1} = eqPoint1;
eqPoints{2} = eqPoint2;

%
%
endTime = max(eqPoint1(1,end),eqPoint2(1,end)) + 20;
startTime = 0;
stepForSolve = 0.01;
stepForOutput = 2; %(s)
tspan = startTime:stepForOutput:endTime;
%
%
usePreparedPsai = 0;%%%%% note

Psai = [
%     -33.00 -33.00 -33.00 -33.00 -33.00 -33.00;
    -88.00 -33.00 -33.00 -88.00 -33.00 -33.00;
    -75.00 -75.00 -33.00 -75.00 -75.00 -33.00;
    -30.00 -85.00 -30.00 -30.00 -85.00 -30.00;
    -30.00 -75.00 -75.00 -30.00 -75.00 -75.00;
    -30.00 -30.00 -85.00 -30.00 -30.00 -85.00;
    -75.00 -30.00 -75.00 -75.00 -30.00 -75.00;
    -85.00 -30.00 -30.00 -85.00 -30.00 -30.00;
    -33.00 -33.00 -33.00 -33.00 -33.00 -33.00;
    -33.00 -33.00 -33.00 -33.00 -33.00 -33.00;    
    ]'*(pi/180);


configs.simulations.startTime = startTime;
configs.simulations.endTime = endTime;
configs.simulations.stepForSolve = stepForSolve;
configs.simulations.stepForOutput = stepForOutput;
configs.simulations.tspan = tspan;
configs.eqPoints = eqPoints;
configs.usePreparedPsai = usePreparedPsai;
configs.Psai = Psai;


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

if exist('x_fp_0', 'var') == 1
    configs.fps.x_fp_0 = x_fp_0;
    configs.fps.y_fp_0 = y_fp_0;
    configs.fps.t_fp_0 = t_fp_0;
    configs.fps.m_fp_0 = m_fp_0;
    configs.fps.i_fp_0 = i_fp_0;
    configs.fps.fps = fps;
else
    configs.fps.x_fp_0 = [];
    configs.fps.y_fp_0 = [];
    configs.fps.t_fp_0 = [];
    configs.fps.m_fp_0 = [];
    configs.fps.i_fp_0 = [];
    configs.fps.fps = [];
end


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


if exist('walls', 'var') == 1
    configs.walls = walls;
else
    configs.walls = [];
end


%%

a = stepOne(handles, args, configs, plotOptions);
