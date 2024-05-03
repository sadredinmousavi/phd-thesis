 function [solverInputs, solverOutputs] = solveTheSystem(configs, args, savingData)

%% Cylinder magnetic field visualization
 % This script is a DEMO for the visualization of the magnetic flux density 
 % field lines of an axially magnetized cylinder. 


%% 
% plotOptions = vars.plotOptions.dynamic;

options = odeset('OutputFcn',@odeprog,'Events',@odeabort);
ans0_mr = [configs.mrs.x_mr_0 configs.mrs.y_mr_0 zeros(1,length(configs.mrs.x_mr_0)) zeros(1,length(configs.mrs.y_mr_0))];
ans0_fp = [configs.fps.x_fp_0 configs.fps.y_fp_0 zeros(1,length(configs.fps.x_fp_0)) zeros(1,length(configs.fps.y_fp_0)) configs.fps.t_fp_0 zeros(1,length(configs.fps.t_fp_0))];
ans0    = [ans0_mr ans0_fp];
%
% inputs = vars.dynamicSolverInputs;
inputs.walls  = configs.walls;
inputs.mr_num = length(configs.mrs.x_mr_0);
inputs.fp_num = length(configs.fps.x_fp_0);
inputs.args = args;
inputs.fps  = configs.fps.fps;
inputs.m_fp = configs.fps.m_fp_0;
inputs.i_fp = configs.fps.i_fp_0;
inputs.k_fp = configs.fps.k_fp_0;
inputs.r_mr = configs.mrs.r_mr_0;
inputs.m_mr = configs.mrs.m_mr_0;
inputs.k_mr = configs.mrs.k_mr_0;
inputs.drag_coeff = args.drag_coeff;
inputs.viscosity = args.viscosity;
inputs.sigma = args.sigma;
inputs.epsilun = args.epsilun;
inputs.threshold = args.threshold;
%
% designPsaiController(vars);
[t,ans1] = ode45(@(t,y)systemDynamics(t,y,inputs), configs.simulations.tspan, ans0, options);
% [t,ans1] = myRungeKutta(@(t,y)systemDynamics(t,y,inputs), vars.tspan, ans0);
Npoints = length(configs.vars.x_space);
%
%
waitText  = 'Gathering plot data  - Please wait...';
waitIters = size(ans1,1);
waitHandle = waitbar(0,waitText);
for i=1:size(ans1,1)
    waitbar(i/waitIters,waitHandle, waitText);
    [Psai,eqPoints_sequence,targetInd,real_eqPoint_x_seq,real_eqPoint_y_seq] = psaiController(t(i));
%     Psai = psaiController(t(i));
%     if plotOptions.fieldVectors
        [~, Frho, Faxial] = calculateForceField(configs.vars.x_space,configs.vars.y_space,Psai);
        plotData(i).Frho = Frho;
        plotData(i).Faxial = Faxial;
%     end
%     plotData(i).eqPoints;
    plotData(i).targetInd = targetInd;
    plotDataStatic.real_eqPoint_x_seq(:,i) = real_eqPoint_x_seq;
    plotDataStatic.real_eqPoint_y_seq(:,i) = real_eqPoint_y_seq;
    plotDataStatic.eqPoints_sequence = eqPoints_sequence;
    plotData(i).Psai = Psai;
end
close(waitHandle)

% cd('data')
% % save(savingData.fileName3, 'vars', 'inputs', 'plotData', 'plotDataStatic', 't', 'ans1')
% save(savingData.fileName3, 'inputs', 'plotData', 'plotDataStatic', 't', 'ans1')
% cd('..')

solverInputs = inputs;
solverOutputs.plotData = plotData;
solverOutputs.plotDataStatic = plotDataStatic;
solverOutputs.t = t;
solverOutputs.ans = ans1;



function [tspan, X] = myRungeKutta(Fun, tspan, X0)
    waitText_  = 'Solving the equations  - Please wait...';
    waitIters_ = length(tspan);
    waitHandle_ = waitbar(0,waitText_);
    dt = tspan(2)-tspan(1);
    X(:,1)  = X0;
    for i_=1:length(tspan)-1
        waitbar(i_/waitIters_,waitHandle_, waitText_);
        t_  = tspan(i_);
        Xi = X(:,i_);
        K1 = Fun(t_,Xi);
        K2 = Fun(t_+dt/2,Xi+K1*dt/2);
        K3 = Fun(t_+dt/2,Xi+K2*dt/2);
        K4 = Fun(t_+dt,Xi+K3*dt);
        Xi = Xi+(K1+2*K2+2*K3+K4)/6*dt;
        X(:,i_+1)=Xi;
    end
    X = X';
    close(waitHandle_)
end


end

