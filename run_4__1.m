%% Cylinder magnetic field visualization
 % This script is a DEMO for the visualization of the magnetic flux density 
 % field lines of an axially magnetized cylinder. 


%%
plotOptions = vars.plotOptions.dynamic;

options=odeset('OutputFcn',@odeprog,'Events',@odeabort);
ans0 = [vars.x_mr_0 vars.y_mr_0 zeros(1,length(vars.x_mr_0)) zeros(1,length(vars.x_mr_0))];
[t,ans1] = ode45(@systemDynamics, vars.tspan, ans0, options);
Npoints = length(vars.x_space);
%
%
waitText  = 'Gathering plot data  - Please wait...';
waitIters = size(ans1,1);
waitHandle = waitbar(0,waitText);
for i=1:size(ans1,1)
    waitbar(i/waitIters,waitHandle, waitText);
    [Psai,eqPoint1,eqPoint2,targetInd,real_eqPoint_x_seq,real_eqPoint_y_seq] = psaiController(t(i));
%     Psai = psaiController(t(i));
    if plotOptions.fieldVectors
        [~, Frho, Faxial] = calculateForceField(vars.x_space,vars.y_space,Psai);
        plotData(i).Frho = Frho;
        plotData(i).Faxial = Faxial;
    end
    plotData(i).eqPoint1 = eqPoint1;
    plotData(i).eqPoint2 = eqPoint2;
    plotData(i).targetInd = targetInd;
    plotDataStatic.real_eqPoint_x_seq(:,i) = real_eqPoint_x_seq;
    plotDataStatic.real_eqPoint_y_seq(:,i) = real_eqPoint_y_seq;
    plotData(i).Psai = Psai;
end
close(waitHandle)



%%
p1 = figure;
set(p1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf, 'color', 'w');
% pause(2)
cd('data')




for i=1:size(ans1,1)
    n = size(ans1,2)/4;
    delete(findall(gcf,'type','annotation'))
    plot(ans1(i,1:n),ans1(i,n+1:2*n),'bo', 'MarkerSize',3 , 'LineWidth',3 );%plot(ans1(i,1:n),ans1(i,n+1:2*n),'bo', 'MarkerSize',8 , 'LineWidth',3 );
    hold on
    if plotOptions.robotsTrack
        for j=1:n
            plot(ans1(1:i,j),ans1(1:i,n+j),'r-', 'LineWidth',1 );% plot(ans1(1:i,j),ans1(1:i,n+j),'r-', 'LineWidth',3 );
        end
    end
    xlim([-vars.plotDomain vars.plotDomain]);
    ylim([-vars.plotDomain vars.plotDomain]);
    axis square
    %
    %
    if plotOptions.fieldVectors
        plot_field = streamslice(vars.x_space,vars.y_space,plotData(i).Frho,plotData(i).Faxial,'method','cubic');
        set(plot_field,'Color','black','LineWidth',1.2);
    end
    %
    if plotOptions.magnets
        plot(vars.MagPos(:,1), vars.MagPos(:,2), 'g.', 'MarkerSize', 15);
    end
    %
    %
    if plotOptions.areaBorders
        plot([min(x_space) max(x_space) max(x_space) min(x_space) min(x_space)], [min(y_space) min(y_space) max(y_space) max(y_space) min(y_space)], 'k-')
    end
    %
    if plotOptions.eqPoints %real eqpoints is the target when vars.findEqFromMinimization = 0;
        plot(plotDataStatic.real_eqPoint_x_seq(:,i), plotDataStatic.real_eqPoint_y_seq(:,i), 'rx', 'MarkerSize', 25, 'LineWidth', 2);
    end
    if plotOptions.eqPointsTrack
        for cnt=1:size(plotDataStatic.real_eqPoint_x_seq,1)
            plot(plotDataStatic.real_eqPoint_x_seq(cnt,1:i), plotDataStatic.real_eqPoint_y_seq(cnt,1:i), 'r-', 'LineWidth',1 );        
        end
%         plot(plotData(i).target_x_seq(plotData(i).targetInd), plotData(i).target_y_seq(plotData(i).targetInd), 'b.', 'MarkerSize', 25)
    end
    %
    if plotOptions.printPsaiValues + plotOptions.printLambdaValues > 1
        text = {};
        if plotOptions.printPsaiValues
            for k=1:length(plotData(i).Psai)
                text{k} = sprintf('{\\psi}_%d=%.0f^{\\circ}\n', k, plotData(i).Psai(k)*(180/pi));
            end
        end
        if plotOptions.printLambdaValues
            cnt = length(text);
%             for k=(length(text)+1):(vars.eqPointsNum+length(text)+1)
                eqPoint = plotData(i).eqPoint1;
                [r1, a1, b1, c1] = calculateParamsFromPoint(eqPoint, vars.MagPos);
                [isStable1,hessian.point1] = isHessianStable(Psai, a1, b1, c1);
                [V,D] = eig(hessian.point1); % V(:,i)
                angle1 = atan(V(2,1)/V(1,1))*(180/pi);
                angle2 = atan(V(2,2)/V(1,2))*(180/pi);
                d1 = sqrt(abs(D(1,1)));
                d2 = sqrt(abs(D(2,2)));
                k=1;
                text{cnt+1} = sprintf('point%d:[{\\theta}_1=%.0f^{\\circ}, {\\theta}_2=%.0f^{\\circ}]\n', k, angle1, angle2);
                text{cnt+2} = sprintf('point%d:[{\\lambda}_1=%.0f, {\\lambda}_2=%.0f]\n', k, d1, d2);
                cnt = cnt + 2;
                %
                eqPoint = plotData(i).eqPoint2;
                [r1, a1, b1, c1] = calculateParamsFromPoint(eqPoint, vars.MagPos);
                [isStable1,hessian.point1] = isHessianStable(Psai, a1, b1, c1);
                [V,D] = eig(hessian.point1); % V(:,i)
                angle1 = atan(V(2,1)/V(1,1))*(180/pi);
                angle2 = atan(V(2,2)/V(1,2))*(180/pi);
                d1 = sqrt(abs(D(1,1)));
                d2 = sqrt(abs(D(2,2)));
                k=2;
                text{cnt+1} = sprintf('point%d:[{\\theta}_1=%.0f^{\\circ}, {\\theta}_2=%.0f^{\\circ}]\n', k, angle1, angle2);
                text{cnt+2} = sprintf('point%d:[{\\lambda}_1=%.0f, {\\lambda}_2=%.0f]\n', k, d1, d2);
                cnt = cnt + 2;
%             end
        end
        annotation( 'textbox', 'String', sprintf('%s', text{:}), 'Color', 'k', ...
            'FontSize', 14, 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.75,0.9,0.5,0] )
    end
    %
    title(sprintf('t = %1.3f s', t(i)),'Color','r', 'FontSize', 30)
    %
    hold off
    xlabel('x [m]','interpreter','latex')
    ylabel('y [m]','interpreter','latex')
    set(gca, 'fontsize', 30)
    drawnow
    frame = getframe(p1);
    im{i} = frame2im(frame);
    [A,map] = rgb2ind(im{i}, 256);
    if i==1
        imwrite(A,map,vars.fileName2, 'gif', 'LoopCount', Inf, 'DelayTime', t(2)-t(1));
    else
        imwrite(A,map,vars.fileName2, 'gif', 'WriteMode', 'append', 'DelayTime', t(2)-t(1));
    end
end
cd('..')

