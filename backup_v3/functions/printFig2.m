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
    
  
    plot_field = streamslice(vars.x_space,vars.y_space,Frho,Faxial,'method','cubic');
    set(plot_field,'Color','black','LineWidth',1.2);
    hold on 
    plot(vars.MagPos(:,1), vars.MagPos(:,2), 'r.', 'MarkerSize', 15)
    %
    if options.plotEqPoints
        for j=1:length(options.eq_points)
            plot(options.eq_points{j}(1), options.eq_points{j}(2), 'rx', 'MarkerSize', 15,'LineWidth',2);
        end
    end
    if options.printPsaiValues + options.printLambdaValues > 1
        text = {};
        if options.printPsaiValues
            for k=1:length(Psai)
                text{k} = sprintf('{\\psi}_%d=%.0f^{\\circ}\n', k, Psai(k)*(180/pi));
            end
        end
        if options.printLambdaValues
            cnt = length(text);
            for k=1:length(options.eq_points)
                eqPoint = options.eq_points{k};
                [r1, a1, b1, c1] = calculateParamsFromPoint(eqPoint, vars.MagPos);
                [isStable1,hessian.point1] = isHessianStable(Psai, a1, b1, c1);
                [V,D] = eig(hessian.point1); % V(:,i)
                angle1 = atan(V(2,1)/V(1,1))*(180/pi);
                angle2 = atan(V(2,2)/V(1,2))*(180/pi);
                d1 = sqrt(abs(D(1,1)));
                d2 = sqrt(abs(D(2,2)));
                text{cnt+1} = sprintf('point%d:[{\\theta}_1=%.0f^{\\circ}, {\\theta}_2=%.0f^{\\circ}]\n', k, angle1, angle2);
                text{cnt+2} = sprintf('point%d:[{\\lambda}_1=%.0f, {\\lambda}_2=%.0f]\n', k, d1, d2);
                cnt = cnt + 2;
            end
        end
        annotation( 'textbox', 'String', sprintf('%s', text{:}), 'Color', 'k', ...
            'FontSize', 14, 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.75,0.9,0.5,0] )
    end
    %
    title('\textbf{Force Field}','interpreter','latex', 'fontsize',14)
    xlabel('x [m]','interpreter','latex')
    ylabel('y [m]','interpreter','latex')
    xlim([-vars.plotDomain vars.plotDomain]);
    ylim([-vars.plotDomain vars.plotDomain]);
    axis square

    hold off
    
    if strcmp(saveFigName, '') == 0
        cd('print')
        set(p1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % saveas(p1, '1aa.png')
        print(p1,saveFigName,'-dpng','-r300')
        close all
        cd('..')
    end
    
    
    
    

end

