d = stepOne("10_object.mat");
% pso = [0; 90; 90; 90; 90; 90]*(pi/180);
% force_field_symbolic(0,0,pso)

% d.args.pm.L = 40 * 1e-3; %[m]
% d.args.pm.D = 40 * 1e-3; %[m]
% d.args.pm.V = d.args.pm.L ^ 3;
% d.args.pm.m = d.args.M * d.args.pm.V;
% d.symbolicCalculations();
% d.saveToFile();

% d.handles.calcPsaiFromEqFunc = @calculatePsai_minimization;%@calculatePsai_6PM; %@calculatePsai_4PM; %@calculatePsai_minimization;
% d.handles.findEqFromMinimization = 0;
% d.saveToFile
% d.designPsai

% Psai = [-63;-33;-33;-33;-33;-33] * (pi/180);
% Psai = d.designPsaiOutput.Psai_t(:,1);
% d.calcEqPointByDynamics(Psai)


eqPoints{1} = convertEqPointRelativeTimeToAbs( [ [0;+0.02;+0.04;15] [30;-0.02;+0.04;-1] [10;-0.02;+0.04;-1] [10;-0.00;+0.06;-1] [10;+0.00;+0.04;-1] [10;+0.00;+0.02;-1] [10;+0.00;-0.06;-1] [10;-0.02;-0.06;-1] [10;-0.06;-0.06;-1] ]);
eqPoints{2} = convertEqPointRelativeTimeToAbs( [ [0;+0.02;+0.04;15] [30;+0.02;+0.04;-1] [10;-0.02;+0.04;-1] [10;-0.00;+0.06;-1] [10;+0.00;+0.04;-1] [10;+0.00;+0.02;-1] [10;+0.00;-0.06;-1] [10;+0.02;-0.06;-1] [10;+0.06;-0.06;-1] ]);
d.setEqPoints(eqPoints);
d.designPsai();
d.saveToFile();


% walls_1 = [ [0;-0.04;+0.04;-0.03;-0.04] [0;-0.03;-0.04;-0.03;-0.09] [0;-0.03;-0.09;+0.09;-0.09] ];
% % walls_2 = [ [0;+0.04;+0.04;+0.03;-0.02] [0;+0.03;-0.02;+0.04;-0.03] [0;+0.04;-0.03;+0.09;-0.03] ];
% walls_2 = [ [0;+0.04;+0.04;+0.03;-0.04] [0;+0.03;-0.04;+0.09;-0.04] ];
% % walls_1 = [ [0;-0.06;+0.08;-0.02;+0.04] [0;-0.02;+0.04;+0.02;-0.00] [0;+0.02;-0.00;+0.06;-0.04] ];
% % walls_2 = [ [0;+0.06;-0.04;+0.06;+0.04] ];
% walls = [ walls_1 walls_2];
% d.setWalls(walls);
% d.saveToFile();

clear seedData
seedData{1}.center = [-0.06, +0.10];
seedData{1}.num = 5; %4
seedData{1}.offset = 0.02;
%
seedData{2}.center = [+0.06, +0.10];
seedData{2}.num = 5; %4
seedData{2}.offset = 0.02;
%
d.setMRsLoc(seedData)
d.saveToFile();


clear fps
fps = {};
%
d.setFPsLoc(fps)
d.saveToFile();






% d.symbolicCalculations();
% d.designPsai();
d.runSimulations();
d.saveToFile();

plotOptionsDyn.fieldVectors = 1;
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