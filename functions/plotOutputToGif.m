function [] = plotOutputToGif(solverInputs, solverOutputs, configs, savingData, plotOptions)




p1 = figure;
set(p1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf, 'color', 'w');
% pause(2)
% cd('data')
%
% load(dataFileName)
% plotOptions = vars.plotOptions.dynamic;
% if nargin > 1
%     plotOptions = plotOptionsFromInput;
% end
%
n = solverInputs.mr_num;
m = solverInputs.fp_num;
w = size(solverInputs.walls,2);
ans1           = solverOutputs.ans;
t              = solverOutputs.t;
plotData       = solverOutputs.plotData;
plotDataStatic = solverOutputs.plotDataStatic;
for i=1:size(ans1,1)
    delete(findall(gcf,'type','annotation'))
    clf(gcf)
    title(sprintf('t = %1.3f s', t(i)),'Color','r', 'FontSize', 30)
    hold on
    % microrobots
%     plot(ans1(i,1:n),ans1(i,n+1:2*n),'bo', 'MarkerSize',3 , 'LineWidth',3 );
    for j=1:n
        x_mr = ans1(i,j);
        y_mr = ans1(i,n+j);
        r_mr = solverInputs.r_mr(j)*3; %%%note
        position = [x_mr-r_mr, y_mr-r_mr, 2*r_mr, 2*r_mr];
        rectangle('Position',position,'Curvature',1,'FaceColor','b','EdgeColor','none')
    end
    if plotOptions.robotsTrack
        for j=1:n
            plot(ans1(1:i,j),ans1(1:i,n+j),'r-', 'LineWidth',1 );
        end
    end 
    % particles
    if m > 0
        for j=1:m
            fp = solverInputs.fps{j};
            x_fp = ans1(i,4*n+j);
            y_fp = ans1(i,4*n+m+j);
            if fp.type == 1
                r_fp = fp.radius;
                position = [x_fp-r_fp, y_fp-r_fp, 2*r_fp, 2*r_fp];
                rectangle('Position',position,'Curvature',1,'FaceColor','r','EdgeColor','none')
            elseif fp.type == 2
                t_fp = ans1(i,4*n+4*m+j);
                fp_length = fp.length;
                position = [x_fp-fp_length/2, y_fp-fp_length/2, fp_length, fp_length];
                rectangle('Position',position,'Curvature',0,'FaceColor','r','EdgeColor','none')
            elseif fp.type == 3
                t_fp = ans1(i,4*n+4*m+j);
                fp_points = fp.points + [x_fp; y_fp];
                plot([fp_points(1,:) fp_points(1,1)], [fp_points(2,:) fp_points(2,1)],'r-', 'LineWidth',1 );
            end
            if plotOptions.particlesTrack
                plot(ans1(1:i,4*n+j),ans1(1:i,4*n+m+j),'r-', 'LineWidth',1 );
            end
        end
    end
    %
    if plotOptions.fieldVectors
        plot_field = streamslice(configs.vars.x_space,configs.vars.y_space,plotData(i).Frho,plotData(i).Faxial,'method','cubic');
        set(plot_field,'Color','black','LineWidth',1.2);
    end
    %
    if plotOptions.magnets
        plot(configs.epms.MagPos(:,1), configs.epms.MagPos(:,2), 'g.', 'MarkerSize', 15);
    end
    %
    if plotOptions.channelLines
        if w > 0
            for cnt=1:w
                plot(solverInputs.walls([2 4], cnt), solverInputs.walls([3 5], cnt),'k-');
            end
        end
    end
    %
    if plotOptions.areaBorders
        plot([min(configs.vars.x_space) max(configs.vars.x_space) max(configs.vars.x_space) min(configs.vars.x_space) min(configs.vars.x_space)], [min(configs.vars.y_space) min(configs.vars.y_space) max(configs.vars.y_space) max(configs.vars.y_space) min(configs.vars.y_space)], 'k-')
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
            for k=1:size(plotDataStatic.real_eqPoint_x_seq,1)
                eqPoint = [plotDataStatic.real_eqPoint_x_seq(k,i) plotDataStatic.real_eqPoint_y_seq(k,i)];
                [r1, a1, b1, c1] = calculateParamsFromPoint(eqPoint, configs.epms.MagPos);
                [isStable1,hessian.point1] = isHessianStable(plotData(i).Psai, a1, b1, c1);
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
    hold off
    xlim([-configs.vars.plotDomain configs.vars.plotDomain]);
    ylim([-configs.vars.plotDomain configs.vars.plotDomain]);
    axis square
    xlabel('x [m]','interpreter','latex')
    ylabel('y [m]','interpreter','latex')
    set(gca, 'fontsize', 30)
    drawnow
    frame = getframe(p1);
    im{i} = frame2im(frame);
    [A,map] = rgb2ind(im{i}, 256);
    if i==1
        imwrite(A,map,savingData.fileName2, 'gif', 'LoopCount', Inf, 'DelayTime', t(2)-t(1));
    else
        imwrite(A,map,savingData.fileName2, 'gif', 'WriteMode', 'append', 'DelayTime', t(2)-t(1));
    end
end
% cd('..')


end

