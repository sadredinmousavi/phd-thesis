function [] = solveTheSystem(vars)

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
[t,ans1] = ode15s(@(t,y)systemDynamics(t,y,inputs), vars.tspan, ans0, options);
% [t,ans1] = myRungeKutta(@(t,y)systemDynamics(t,y,inputs), vars.tspan, ans0);
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
%     if plotOptions.fieldVectors
        [~, Frho, Faxial] = calculateForceField(vars.x_space,vars.y_space,Psai);
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

cd('data')
save(vars.fileName3, 'vars', 'inputs', 'plotData', 'plotDataStatic', 't', 'ans1')
cd('..')


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

