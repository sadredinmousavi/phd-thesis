%% Cylinder magnetic field visualization
 % This script is a DEMO for the visualization of the magnetic flux density 
 % field lines of an axially magnetized cylinder. 


%% 
plotOptions = vars.plotOptions.dynamic;

options = odeset('OutputFcn',@odeprog,'Events',@odeabort);
ans0_mr = [vars.x_mr_0 vars.y_mr_0 zeros(1,length(vars.x_mr_0)) zeros(1,length(vars.y_mr_0))];
ans0_fp = [vars.x_fp_0 vars.y_fp_0 zeros(1,length(vars.x_fp_0)) zeros(1,length(vars.y_fp_0)) vars.t_fp_0 zeros(1,length(vars.t_fp_0))];
ans0    = [ans0_mr ans0_fp];
%
% inputs = vars.dynamicSolverInputs;
inputs.walls  = vars.walls;
inputs.mr_num = length(vars.x_mr_0);
inputs.fp_num = length(vars.x_fp_0);
inputs.args = vars.args;
inputs.fps = vars.fps;
inputs.m_fp = vars.m_fp_0;
inputs.i_fp = vars.i_fp_0;
inputs.r_mr = vars.r_mr_0;
inputs.m_mr = vars.m_mr_0;
inputs.drag_coeff = vars.drag_coeff;
inputs.sigma = vars.sigma;
inputs.epsilun = vars.epsilun;
inputs.threshold = vars.threshold;
%
designPsaiController(vars);
[t,ans1] = ode45(@(t,y)systemDynamics(t,y,inputs), vars.tspan, ans0, options);
Npoints = length(vars.x_space);
%
%
waitText  = 'Gathering plot data  - Please wait...';
waitIters = size(ans1,1);
waitHandle = waitbar(0,waitText);
for i=1:size(ans1,1)
    waitbar(i/waitIters,waitHandle, waitText);
    [Psai,eqPoints_sequence,targetInd,real_eqPoint_x_seq,real_eqPoint_y_seq] = psaiController(t(i));
%     Psai = psaiController(t(i));
    if plotOptions.fieldVectors
        [~, Frho, Faxial] = calculateForceField(vars.x_space,vars.y_space,Psai);
        plotData(i).Frho = Frho;
        plotData(i).Faxial = Faxial;
    end
%     plotData(i).eqPoints;
    plotData(i).targetInd = targetInd;
    plotDataStatic.real_eqPoint_x_seq(:,i) = real_eqPoint_x_seq;
    plotDataStatic.real_eqPoint_y_seq(:,i) = real_eqPoint_y_seq;
    plotDataStatic.eqPoints_sequence = eqPoints_sequence;
    plotData(i).Psai = Psai;
end
close(waitHandle)

cd('data')
save(vars.fileName3, 'vars', 'inputs', 'plotData', 'plotDataStatic', 't', 'ans1')
cd('..')

plotOutputToGif(vars.fileName3)

%%


