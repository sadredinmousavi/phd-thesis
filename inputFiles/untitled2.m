d = stepOne("1_object.mat");
% d.handles.findEqFromMinimization = 0;
% d.saveToFile

% d.designPsai
% Psai = [-63;-33;-33;-33;-33;-33] * (pi/180);
% d.calcEqPointByDynamics(Psai)

% d.runSimulations





% eqPoints{1} = [ [0;-0.02;+0.08] [10;-0.00;+0.06] [20;-0.00;+0.04] [30;-0.00;+0.00] [40;+0.00;-0.04] [50;+0.04;-0.06] [60;+0.08;-0.06] ];
% eqPoints{2} = [ [0;+0.02;+0.08] [10;+0.00;+0.06] [20;+0.00;+0.04] [30;+0.00;+0.00] [40;+0.00;-0.04] [50;+0.04;-0.06] [60;+0.08;-0.06]];
% d.setEqPoints(eqPoints);
% d.designPsai();
% d.saveToFile();


% walls_1 = [ [0;-0.04;+0.04;-0.03;-0.04] [0;-0.03;-0.04;-0.03;-0.09] [0;-0.03;-0.09;+0.09;-0.09] ];
% % walls_2 = [ [0;+0.04;+0.04;+0.03;-0.02] [0;+0.03;-0.02;+0.04;-0.03] [0;+0.04;-0.03;+0.09;-0.03] ];
% walls_2 = [ [0;+0.04;+0.04;+0.03;-0.04] [0;+0.03;-0.04;+0.09;-0.04] ];
% % walls_1 = [ [0;-0.06;+0.08;-0.02;+0.04] [0;-0.02;+0.04;+0.02;-0.00] [0;+0.02;-0.00;+0.06;-0.04] ];
% % walls_2 = [ [0;+0.06;-0.04;+0.06;+0.04] ];
% walls = [ walls_1 walls_2];
% d.setWalls(walls);
% d.saveToFile();

% clear seedData
% seedData{1}.center = [-0.03, +0.07];
% seedData{1}.num = 3; %4
% seedData{1}.offset = 0.02;
% %
% seedData{2}.center = [+0.03, +0.07];
% seedData{2}.num = 3; %4
% seedData{2}.offset = 0.02;
% %
% d.setMRsLoc(seedData)
% d.saveToFile();


% clear fps
% fps{1}.type = 1;
% fps{1}.center = [+0.00 +0.06];
% fps{1}.theta0 = +0.00;
% fps{1}.radius = 0.006 ;
% fps{1}.height = 0.002 ;
% fps{1}.K = [16/3 16/3] * r_fp ;
% %
% d.setFPsLoc(fps)
% d.saveToFile();


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