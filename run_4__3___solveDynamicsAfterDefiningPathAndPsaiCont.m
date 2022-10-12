%% Cylinder magnetic field visualization
 % This script is a DEMO for the visualization of the magnetic flux density 
 % field lines of an axially magnetized cylinder. 


%%


ans0 = [vars.x_mr_0 vars.y_mr_0 zeros(1,length(vars.x_mr_0)) zeros(1,length(vars.x_mr_0))];
[t,ans1] = ode45(@systemDynamics, vars.tspan, ans0);
Npoints = length(vars.x_space);
%
%
waitText  = 'plotting force field - Please wait...';
waitIters = size(ans1,1);
waitHandle = waitbar(0,waitText);
for i=1:size(ans1,1)
    waitbar(i/waitIters,waitHandle, waitText);
    [Psai,eqPoint1,eqPoint2,targetInd,real_eqPoint_x_seq,real_eqPoint_y_seq] = psaiController(t(i));
%     Psai = psaiController(t(i));
    [~, Frho, Faxial] = calculateForceField(vars.x_space,vars.y_space,Psai);
    forceField(i).Frho = Frho;
    forceField(i).Faxial = Faxial;
    forceField(i).eqPoint1 = eqPoint1;
    forceField(i).eqPoint2 = eqPoint2;
    forceField(i).targetInd = targetInd;
    forceField(i).real_eqPoint_x_seq = real_eqPoint_x_seq;
    forceField(i).real_eqPoint_y_seq = real_eqPoint_y_seq;
    forceField(i).Psai = Psai;
end
close(waitHandle)



%%
p1 = figure(2);
set(p1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf, 'color', 'w');
% pause(2)
cd('data')


for i=1:size(ans1,1)
    n = size(ans1,2)/4;
    delete(findall(gcf,'type','annotation'))
    plot(ans1(i,1:n),ans1(i,n+1:2*n),'bo', 'MarkerSize',3 , 'LineWidth',3 );%plot(ans1(i,1:n),ans1(i,n+1:2*n),'bo', 'MarkerSize',8 , 'LineWidth',3 );
    hold on
    for j=1:n
        plot(ans1(1:i,j),ans1(1:i,n+j),'r-', 'LineWidth',1 );% plot(ans1(1:i,j),ans1(1:i,n+j),'r-', 'LineWidth',3 );
    end
    xlim([-vars.plotDomain vars.plotDomain]);
    ylim([-vars.plotDomain vars.plotDomain]);
    axis square
    %
    %
%     hold on
    plot_field = streamslice(vars.x_space,vars.y_space,forceField(i).Frho,forceField(i).Faxial,'method','cubic');
    set(plot_field,'Color','black','LineWidth',1.2);
    %
    plot(vars.MagPos(:,1), vars.MagPos(:,2), 'g.', 'MarkerSize', 15)
    %
%     plot(forceField(i).target_x_seq, forceField(i).target_y_seq, 'b-', 'LineWidth',3 );
%     plot(forceField(i).real_eqPoint_x_seq, forceField(i).real_eqPoint_y_seq, 'b.', 'MarkerSize', 25)
%     plot(forceField(i).target_x_seq(forceField(i).targetInd), forceField(i).target_y_seq(forceField(i).targetInd), 'b.', 'MarkerSize', 25)
    %
    text = {};
    for k=1:length(forceField(i).Psai)
        text{k} = sprintf('{\\psi}_%d=%.0f^{\\circ}\n', k, forceField(i).Psai(k)*(180/pi));
    end
    annotation( 'textbox', 'String', sprintf('%s', text{:}), 'Color', 'k', ...
        'FontSize', 14, 'Units', 'normalized', 'EdgeColor', 'none', ...
        'Position', [0.75,0.9,0.5,0] )
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

