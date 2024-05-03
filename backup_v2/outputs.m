% plotOutputToGif('10_datas.mat');%ode45 free
% plotOutputToGif('22_datas.mat');%ode113
% plotOutputToGif('23_datas.mat');%ode15s
% plotOutputToGif('24_datas.mat');%ode23
% plotOutputToGif('25_datas.mat');%ode23s
% plotOutputToGif('27_datas.mat');%ode45
%
% plotOutputToGif('29_datas.mat');%ode15s 100*sigma
% plotOutputToGif('30_datas.mat');%ode15  100*sigma
% plotOutputToGif('31_datas.mat');%ode15  100*sigma
% plotOutputToGif('31_datas.mat');%ode15  100*sigma
%
% plotOutputToGif('41_datas.mat');%ode45  100*sigma without fp old path
% plotOutputToGif('43_datas.mat');%ode45  100*sigma without fp new path
% plotOutputToGif('44_datas.mat');%ode45  100*sigma without fp new path 2
% plotOutputToGif('47_datas.mat');%ode45  100*sigma 6mm circle
% plotOutputToGif('58_datas.mat');%ode45  100*sigma 4mm rectangle

% plotOutputToGif('71_datas.mat');%ode45  100*sigma 6mm circle
% plotOutputToGif('74_datas.mat');%ode45  5*sigma 6mm circle







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

plotOutputToGif('6_datas.mat', plotOptionsDyn);%42% 43
% plotOutputToGif(vars.fileName3, plotOptionsDyn);