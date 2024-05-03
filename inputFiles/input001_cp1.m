% d = stepOne("8_object.mat");
d = stepOne("111_object_2.mat");

% pso = [0; 90; 90; 90; 90; 90]*(pi/180);
% force_field_symbolic(0,0,pso)

% d.args.pm.L = 40 * 1e-3; %[m]
% d.args.pm.D = 40 * 1e-3; %[m]
% d.args.pm.V = d.args.pm.L ^ 3;
% d.args.pm.m = d.args.M * d.args.pm.V;
% d.symbolicCalculations();
d.args.threshold = 5*d.args.sigma + 0.002; %5*sigma + 0.004;
% d.saveToFile();

% d.handles.findEqFromMinimization = 0;
% d.saveToFile
% d.designPsai

% Psai = [-63;-33;-33;-33;-33;-33] * (pi/180);
% Psai = d.designPsaiOutput.Psai_t(:,3);
% d.calcEqPointByDynamics(Psai)



eqP = addEqPointToSerieByDelta(+0.00, +0.10, 5);
eqP = addEqPointToSerieByDelta(+0.00, -0.01, 5, eqP);
eqP = addEqPointToSerieByDelta(+0.00, -0.01, 10, eqP);
eqP = addEqPointToSerieByDelta(+0.00, -0.01, 10, eqP);
eqP = addEqPointToSerieByDelta(+0.00, -0.01, 10, eqP);
eqP = addEqPointToSerieByDelta(+0.00, -0.01, 10, eqP);
eqP = addEqPointToSerieByDelta(+0.00, -0.01, 10, eqP);
eqP = addEqPointToSerieByDelta(+0.00, -0.01, 10, eqP);
eqP = addEqPointToSerieByDelta(+0.00, -0.01, 10, eqP);
eqP = addEqPointToSerieByDelta(+0.01, -0.00, 15, eqP);
eqP = addEqPointToSerieByDelta(+0.01, -0.00, 10, eqP);
eqP = addEqPointToSerieByDelta(+0.01, -0.00, 10, eqP);
eqP = addEqPointToSerieByDelta(+0.00, +0.01, 15, eqP);
eqP = addEqPointToSerieByDelta(+0.00, +0.01, 10, eqP);
eqP = addEqPointToSerieByDelta(+0.00, +0.01, 10, eqP);
eqP = addEqPointToSerieByDelta(-0.01, -0.01, 10, eqP);
eqP = addEqPointToSerieByDelta(-0.01, -0.01, 10, eqP);
eqP = addEqPointToSerieByDelta(-0.01, -0.01, 10, eqP);
eqP = addEqPointToSerieByDelta(-0.01, -0.01, 10, eqP);
eqP = addEqPointToSerieByDelta(-0.01, -0.01, 10, eqP);
eqP = addEqPointToSerieByDelta(-0.01, -0.01, 10, eqP);
eqP = addEqPointToSerieByDelta(-0.01, -0.01, 10, eqP);
eqP = addEqPointToSerieByDelta(-0.01, -0.01, 10, eqP);
eqP = addEqPointToSerieByDelta(-0.01, -0.01, 10, eqP);
eqP = addEqPointToSerieByDelta(-0.00, -0.00, 40, eqP);
eqPoints{1} = [ [0;-0.00;+0.11] eqP ];
eqPoints{2} = [ [0;+0.00;+0.11] eqP ];
d.setEqPoints(eqPoints);
d.designPsai();
% d.saveToFile();


% walls_1 = [ [0;-0.04;+0.04;-0.03;-0.04] [0;-0.03;-0.04;-0.03;-0.09] [0;-0.03;-0.09;+0.09;-0.09] ];
% % walls_2 = [ [0;+0.04;+0.04;+0.03;-0.02] [0;+0.03;-0.02;+0.04;-0.03] [0;+0.04;-0.03;+0.09;-0.03] ];
% walls_2 = [ [0;+0.04;+0.04;+0.03;-0.04] [0;+0.03;-0.04;+0.09;-0.04] ];
% % walls_1 = [ [0;-0.06;+0.08;-0.02;+0.04] [0;-0.02;+0.04;+0.02;-0.00] [0;+0.02;-0.00;+0.06;-0.04] ];
% % walls_2 = [ [0;+0.06;-0.04;+0.06;+0.04] ];
% walls = [ walls_1 walls_2];
% d.setWalls(walls);
% d.saveToFile();

clear seedData
seedData{1}.center = [-0.06, +0.12];
seedData{1}.num = 5; %4
seedData{1}.offset = 0.02;
%
seedData{2}.center = [+0.06, +0.12];
seedData{2}.num = 5; %4
seedData{2}.offset = 0.02;
%
d.setMRsLoc(seedData)
% d.saveToFile();


clear fps
fps{1}.type = 1;
fps{1}.center = [+0.00 +0.08];
fps{1}.theta0 = +0.00;
fps{1}.radius = 0.014 ;
fps{1}.height = 0.002 ;
fps{1}.K = [16/3 16/3] * d.args.fp.D/2 ;
%
d.setFPsLoc(fps)
% d.saveToFile();





% d.symbolicCalculations();
% d.designPsai();
% d.runSimulations();
% d.saveToFile("111_object_2.mat");

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