function [] = printFig(vars, Psai, hasFigure, saveFigName)

    if nargin < 4
        saveFigName = '';
    end
    L = vars.args.pm.L;
    D = vars.args.pm.D;
    %
    for i=1:length(Psai)
        titleText{i} = sprintf('${\\psi}_%d=%.0f^{\\circ}$  ', i, Psai(i)*(180/pi));
    end
    %
    [~, Frho, Faxial] = calculateForceField(vars.x_space,vars.y_space,Psai);
    [~, Brho, Baxial] = calculateMagneticField(vars.x_space,vars.y_space,Psai);
    options = vars.plotOptions.static;
    
    %%
    if nargin < 3
        p1 = figure;
    else
        if hasFigure ~= 0
            p1 = figure(hasFigure);
        else
            p1 = figure;
        end
    end
    
    subplot(3,2,[1 3]) %%
    hold on
    plot_field = surface(vars.x_space,vars.y_space,sqrt(Brho.^2+Baxial.^2));
%     plot_field = streamslice(vars.x_space,vars.y_space,Brho,Baxial,'method','cubic');
%     set(plot_field,'Color','black','LineWidth',1.2);
%     plot(vars.MagPos(:,1), vars.MagPos(:,2), 'r.', 'MarkerSize', 15)
%     for i = 1:size(vars.MagPos,1)
%         drawCylindricalMagnet(L,D,vars.MagPos(i,:),'texture','axial')
%     end
    title('\textbf{Magnetic Field}','interpreter','latex', 'fontsize',14)
    xlabel('x [m]','interpreter','latex')
    ylabel('y [m]','interpreter','latex')
    xlim([-vars.plotDomain vars.plotDomain]);
    ylim([-vars.plotDomain vars.plotDomain]);
    axis square
    %
    subplot(3,2,[2 4]) %%
    plot_field = streamslice(vars.x_space,vars.y_space,Frho,Faxial,'method','cubic');
    set(plot_field,'Color','black','LineWidth',1.2);
    hold on 
    plot(vars.MagPos(:,1), vars.MagPos(:,2), 'r.', 'MarkerSize', 15)
    for i = 1:size(vars.MagPos,1)
        drawCylindricalMagnet(L,D,vars.MagPos(i,:),'texture','axial')
    end
    %
    if options.plotEqPoints
        plot(options.eq_point1(1), options.eq_point1(2), 'rx', 'MarkerSize', 15,'LineWidth',2);
        plot(options.eq_point2(1), options.eq_point2(2), 'rx', 'MarkerSize', 15,'LineWidth',2);
    end
    %
    title('\textbf{Force Field}','interpreter','latex', 'fontsize',14)
    xlabel('x [m]','interpreter','latex')
    ylabel('y [m]','interpreter','latex')
    xlim([-vars.plotDomain vars.plotDomain]);
    ylim([-vars.plotDomain vars.plotDomain]);
    axis square

    subplot(3,2,[5 6])
    for i = 1:size(vars.MagPos,1)
        hold on
        drawCylindricalMagnet(L,D,vars.MagPos(i,:),'texture','axial')
        hold on
        drawVec(vars.MagPos(i,1:3),vars.MagPos(i,4:6)*0.005)
    end
    title(sprintf('%s', titleText{:}) ,'interpreter','latex', 'fontsize',14)
    hold on
    showFrame([0 0 0],[1 0 0]/100,[0 1 0]/100,[0 0 1]/100,...
                    'arrowheadlength',0.35,...
                    'arrowheadwidth',0.2)
    % [caz,cel] = view
    view(-12,24);
    hold off
    % axis([-80 80 -80 80 -80 80]/1000)
    setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera')
    axis equal

    if strcmp(saveFigName, '') == 0
        cd('print')
        set(p1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % saveas(p1, '1aa.png')
        print(p1,saveFigName,'-dpng','-r300')
        close all
        cd('..')
    end
    
    
    
    

end

