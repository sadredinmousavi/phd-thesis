d = stepOne("15_object");
% pso = [0; 90; 90; 90; 90; 90]*(pi/180);
% force_field_symbolic(0,0,pso)

% d.setEpmsLoc(8, 0.25, 0);
% d.saveToFile();
% 
% d.args.pm.L = 60 * 1e-3; %[m]
% d.args.pm.D = 60 * 1e-3; %[m]
% d.args.pm.V = d.args.pm.L ^ 3;
% d.args.pm.m = d.args.M * d.args.pm.V;
% d.handles.symbolicFunctionHandle = @symbolic_fun_b;
% d.symbolicCalculations();
% d.saveToFile();


% clear eqPoints
% eqP = addEqPointToSerieByDelta(+0.00, -0.08, 10);
% eqP = addEqPointToSerieByDelta(+0.00, -0.01, 10, eqP);
% eqP = addEqPointToSerieByDelta(+0.00, -0.05, 10, eqP);
% eqP = addEqPointToSerieByDelta(+0.00, -0.12, 10, eqP);
% eqPoints{1} = [ eqP ];
% eqPoints{2} = [ eqP ];
clear eqPoints
eqPoints{1} = convertEqPointRelativeTimeToAbs( [ [00;+0.00;+0.00] [30;+0.06;+0.09] [30;+0.06;+0.09] [30;+0.10;+0.09] [30;+0.09;+0.09] [30;+0.14;+0.09] [30;+0.13;+0.04] [30;+0.03;+0.04] [30;-0.06;+0.04] [30;-0.13;+0.04] [30;-0.08;+0.06] [20;-0.04;+0.03] [10;-0.03;+0.04] [10;-0.04;+0.04] [10;-0.05;+0.04] [10;-0.05;+0.03] [15;-0.00;+0.07] [10;-0.01;+0.08] [10;+0.00;+0.08] [10;+0.01;+0.08] [10;+0.01;+0.07] [15;+0.00;+0.00] ]);
eqPoints{2} = convertEqPointRelativeTimeToAbs( [ [00;+0.00;+0.00] [30;+0.05;+0.08] [30;+0.06;+0.09] [30;+0.09;+0.08] [30;+0.09;+0.09] [30;+0.13;+0.08] [30;+0.13;+0.04] [30;+0.03;+0.04] [30;-0.06;+0.04] [30;-0.13;+0.04] [30;-0.08;+0.06] [20;-0.04;+0.03] [10;-0.05;+0.02] [10;-0.04;+0.02] [10;-0.03;+0.02] [10;-0.03;+0.03] [15;-0.00;+0.07] [10;+0.01;+0.06] [10;-0.00;+0.06] [10;-0.01;+0.06] [10;-0.01;+0.07] [15;+0.00;+0.00] ]);
d.setEqPoints(eqPoints);
% d.designPsai();
d.saveToFile();

d.handles.calcPsaiFromEqFunc = @calculatePsai_minimization;%@calculatePsai_minimization_3;%@calculatePsai_6PM; %@calculatePsai_4PM; %@calculatePsai_minimization;
d.handles.findEqFromMinimization = 0;
d.saveToFile
d.designPsai

% Psai = [-63;-33;-33;-33;-33;-33] * (pi/180);
Psai = d.designPsaiOutput.Psai_t(:,7);
d.calcEqPointByDynamics(Psai)


% clear eqPoints
% eqPoints{1} = [ [0;-0.05;+0.08] [50;-0.05;+0.07] [100;-0.03;+0.06] [150;-0.00;+0.05] [200;-0.00;+0.03] [300;+0.00;+0.01] [400;+0.00;-0.01] [500;+0.00;-0.05] [600;+0.00;-0.07] [700;+0.03;-0.09] [800;+0.05;-0.09] [900;+0.09;-0.09] [1000;+0.13;-0.09] [1100;+0.14;-0.06] [1200;+0.14;-0.00] ];
% eqPoints{2} = [ [0;+0.05;+0.08] [50;+0.05;+0.07] [100;+0.03;+0.06] [150;+0.00;+0.05] [200;+0.00;+0.03] [300;+0.00;+0.01] [400;+0.00;-0.01] [500;+0.00;-0.05] [600;+0.00;-0.07] [700;+0.03;-0.09] [800;+0.05;-0.09] [900;+0.09;-0.09] [1000;+0.13;-0.09] [1100;+0.14;-0.06] [1200;+0.14;-0.00] ];
% d.setEqPoints(eqPoints);
% d.designPsai();
% d.saveToFile();


% walls_1 = [ [0;-0.04;+0.04;-0.03;-0.04] [0;-0.03;-0.04;-0.03;-0.09] [0;-0.03;-0.09;+0.09;-0.09] ];
% walls_2 = [ [0;+0.04;+0.04;+0.03;-0.04] [0;+0.03;-0.04;+0.09;-0.04] ];
% walls = [ walls_1 walls_2]; walls = [];
% d.setWalls(walls);
% d.saveToFile();

% clear seedData
% seedData{1}.center = [+0.10, +0.00];
% seedData{1}.num = 6; %4
% seedData{1}.offset = 0.02;
% %
% % seedData{2}.center = [+0.06, +0.10];
% % seedData{2}.num = 5; %4
% % seedData{2}.offset = 0.02;
% %
% d.setMRsLoc(seedData)
% d.saveToFile();


% clear fps
% fps = {};
% %
% d.setFPsLoc(fps)
% d.saveToFile();






% d.symbolicCalculations();
% d.designPsai();
d.runSimulations();
d.saveToFile();

plotOptionsDyn.fieldVectors = 0;
plotOptionsDyn.magnets = 1;
plotOptionsDyn.areaBorders = 1;
plotOptionsDyn.robotsTrack = 0;
plotOptionsDyn.particlesTrack = 0;
plotOptionsDyn.channelLines = 1;
plotOptionsDyn.eqPoints = 1;
plotOptionsDyn.eqPointsTrack = 0;
plotOptionsDyn.printPsaiValues = 1;
plotOptionsDyn.printLambdaValues = 1;

d.plotOutput(plotOptionsDyn)