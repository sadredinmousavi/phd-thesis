% d = stepOne("9_object.mat");
d = stepOne("222_object_1.mat");

% pso = [0; 90; 90; 90; 90; 90]*(pi/180);
% force_field_symbolic(0,0,pso)

% d.setEpmsLoc(8, 0.25, 0);
% d.saveToFile();

% d.args.pm.L = 40 * 1e-3; %[m]
% d.args.pm.D = 40 * 1e-3; %[m]
% d.args.pm.V = d.args.pm.L ^ 3;
% d.args.pm.m = d.args.M * d.args.pm.V;
% d.symbolicCalculations();
% d.saveToFile();

% d.handles.calcPsaiFromEqFunc = @calculatePsai_minimization;%@calculatePsai_6PM; %@calculatePsai_4PM; %@calculatePsai_minimization;
d.handles.findEqFromMinimization = 0;
% d.saveToFile
% d.designPsai




eqPoints{1} = convertEqPointRelativeTimeToAbs( [ [0;-0.05;+0.08] [10;-0.05;+0.07] [10;-0.03;+0.06] [10;-0.00;+0.05] [10;-0.00;+0.03] [10;+0.00;+0.01] [10;+0.00;-0.01] [10;+0.00;-0.05] [10;+0.00;-0.07] [10;+0.03;-0.09] [20;+0.10;-0.07] [35;+0.12;-0.14] ]);
eqPoints{2} = convertEqPointRelativeTimeToAbs( [ [0;+0.05;+0.08] [10;+0.05;+0.07] [10;+0.03;+0.06] [10;+0.00;+0.05] [10;+0.00;+0.03] [10;+0.00;+0.01] [10;+0.00;-0.01] [10;+0.00;-0.05] [10;+0.00;-0.07] [10;+0.03;-0.09] [20;+0.10;-0.07] [35;+0.12;-0.14] ]);
d.setEqPoints(eqPoints);
d.configs.Psai= [
   -1.5231   -1.5255   -1.4990   -0.1765   -0.4890   -0.2287   -0.2433   -0.4649   -0.1246    0.9385      -1.3708  -1.4358
   -1.5482   -1.5540   -1.5177   -1.2520   -0.2473   -0.7581    0.2231   -0.2421    0.0260   -0.1205      -0.4015   0.5326
   -1.4680   -1.4343   -1.1967   -1.1235   -1.4616    0.6344    0.5709   -0.1011    0.0026   -0.0548      -0.3225   1.3158
   -1.5483   -1.5530   -1.5182   -1.2347   -0.2473   -0.7586    0.2228   -0.2446    0.0134   -0.0362       0.1974   0.7453
   -1.5218   -1.5319   -1.4969   -0.2586   -0.4890   -0.2296   -0.2441   -0.4865   -0.1020    0.0503      -0.3225   1.3158
    1.3836    1.2087    0.5488   -0.2536    0.0340    0.2200   -0.7592   -1.2077   -1.4700   -1.5205      -0.4015   0.5326
    0.1777   -1.0542    0.2562   -0.1801    0.1347    0.5675    0.6360   -1.1661   -1.1977   -1.3403      -1.3708  -1.4358
    1.3796    1.2280    0.5402   -0.1538    0.0340    0.2235   -0.7595   -1.2146   -1.4688   -1.5394      -1.4605  -1.5609
];
d.configs.usePreparedPsai = 1;
d.designPsai();
d.saveToFile();


% eqPoints{1} = convertEqPointRelativeTimeToAbs( [ [0;+0.10;+0.00] [0;+0.13;+0.13] [10;+0.07;+0.13]  [25;-0.01;-0.06] [10;+0.08;+0.00] [10;+0.04;+0.00] [10;+0.00;-0.06] [10;-0.02;-0.06] [10;-0.06;-0.06] ]);
% % eqPoints{2} = convertEqPointRelativeTimeToAbs( [ [0;+0.10;+0.00] [0;+0.13;+0.13] [10;+0.07;+0.13]  [25;+0.11;+0.04] [10;+0.08;+0.00] [10;+0.04;+0.00] [10;+0.00;-0.06] [10;-0.02;-0.06] [10;-0.06;-0.06] ]);
% eqPoints{1} = convertEqPointRelativeTimeToAbs( [ [0;-0.05;+0.08] [10;-0.05;+0.07] [10;-0.03;+0.06] [10;-0.00;+0.05] [10;-0.00;+0.03] [10;+0.00;+0.01] [10;+0.00;-0.01] [10;+0.00;-0.05] [10;+0.00;-0.07] [10;+0.03;-0.09] [10;+0.05;-0.09] [10;+0.09;-0.09] [10;+0.10;-0.03] [10;+0.12;-0.00] [10;+0.14;-0.00] ]);
% eqPoints{2} = convertEqPointRelativeTimeToAbs( [ [0;+0.05;+0.08] [10;+0.05;+0.07] [10;+0.03;+0.06] [10;+0.00;+0.05] [10;+0.00;+0.03] [10;+0.00;+0.01] [10;+0.00;-0.01] [10;+0.00;-0.05] [10;+0.00;-0.07] [10;+0.03;-0.09] [10;+0.05;-0.09] [10;+0.09;-0.09] [10;+0.10;-0.03] [10;+0.12;-0.00] [10;+0.14;-0.00] ]);
% % eqPoints{1} = convertEqPointRelativeTimeToAbs( [ [10;-0.00;+0.14] [10;-0.00;+0.10] [10;+0.00;+0.07] [10;+0.00;+0.05] [10;+0.00;+0.03] ]);
% % eqPoints{2} = convertEqPointRelativeTimeToAbs( [ [10;+0.00;+0.14] [10;+0.00;+0.10] [10;+0.00;+0.07] [10;+0.00;+0.05] [10;+0.00;+0.03] ]);
% d.setEqPoints(eqPoints);
% d.designPsai();
% % d.saveToFile();



% % Psai = [-63;-33;-33;-33;-33;-33] * (pi/180);
% for i=1:size(d.designPsaiOutput.Psai_t, 2)
%     Psai = d.designPsaiOutput.Psai_t(:,i)
%     d.calcEqPointByDynamics(Psai)
% %     keyboard;
% end



% walls_1 = [ [0;-0.08;+0.025;-0.02;+0.025] [0;-0.02;+0.025;+0.02;+0.08] ];
% walls_2 = [ [0;-0.08;-0.025;-0.02;-0.025] [0;-0.02;-0.025;+0.02;-0.08] ];
% walls_3 = [ [0;+0.058;+0.08;+0.00;-0.00] [0;+0.00;-0.00;+0.058;-0.08] ];
% % walls_1 = [ [0;-0.06;+0.08;-0.02;+0.04] [0;-0.02;+0.04;+0.02;-0.00] [0;+0.02;-0.00;+0.06;-0.04] ];
% % walls_2 = [ [0;+0.06;-0.04;+0.06;+0.04] ];
% walls = [ walls_1 walls_2 walls_3];
% d.setWalls(walls);
% d.saveToFile();
walls_1 = [ [0;-0.04;+0.04;-0.03;-0.04] [0;-0.03;-0.04;-0.03;-0.09] [0;-0.03;-0.09;+0.09;-0.09] ];
% walls_2 = [ [0;+0.04;+0.04;+0.03;-0.02] [0;+0.03;-0.02;+0.04;-0.03] [0;+0.04;-0.03;+0.09;-0.03] ];
walls_2 = [ [0;+0.04;+0.04;+0.03;-0.04] [0;+0.03;-0.04;+0.09;-0.04] ];
% walls_1 = [ [0;-0.06;+0.08;-0.02;+0.04] [0;-0.02;+0.04;+0.02;-0.00] [0;+0.02;-0.00;+0.06;-0.04] ];
% walls_2 = [ [0;+0.06;-0.04;+0.06;+0.04] ];
walls = [ walls_1 walls_2];
d.setWalls(walls);
d.saveToFile();

clear seedData
seedData{1}.center = [-0.05, +0.08];
seedData{1}.num = 4; %4
seedData{1}.offset = 0.02;
%
seedData{2}.center = [+0.05, +0.08];
seedData{2}.num = 4; %4
seedData{2}.offset = 0.02;
%
d.setMRsLoc(seedData)
d.saveToFile();


% clear fps
% fps = {};
% %
% d.setFPsLoc(fps)
% d.saveToFile();






% d.symbolicCalculations();
% d.designPsai();
d.runSimulations();
d.saveToFile("222_object_1.mat");

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