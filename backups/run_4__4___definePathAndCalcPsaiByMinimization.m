%% Cylinder magnetic field visualization
 % This script is a DEMO for the visualization of the magnetic flux density 

clc, close all, clearvars -except cursor_info
load('input_data')
% load('cursor_info')
%%

if DoHandPreCheck
    
    myPlotter(x,y,Frho,Faxial,MagPos,[], 1, []);
    
    figure
    surf(x,y,F)
      
    return
else
    save('cursor_info','cursor_info')
    cd('data')
    save(fileName4)
    cd('..')
end

%
initialEqPoints = cursor_info(1).Position';
for i=2:length(cursor_info)
    initialEqPoints = [initialEqPoints cursor_info(i).Position'];
end
initialPsai = MagPosR(:,4);
all_eqPoints    = findEqPoints(initialEqPoints,x(2)-x(1),initialPsai);
% myPlotter(x,y,Frho,Faxial,MagPos,initialPsai, 0, eqPoints);
%

%%

eqPointNumber = 1;
dt = 0.2;
endTime = 10;
maxIter = endTime/dt + 1;
maxDisp = 0.1;
target_eqPoint_x(1) = all_eqPoints(1,eqPointNumber);
target_eqPoint_y(1) = all_eqPoints(2,eqPointNumber);
for i=1:maxIter
    target_eqPoint_x(i) = target_eqPoint_x(1) + maxDisp*((i+1)/maxIter);
    target_eqPoint_y(i) = -maxDisp*(((i+1)/maxIter))^4 + target_eqPoint_y(1);
end

%
%

waitText  = 'openloop Psai calculations - Please wait...';
waitHandle = waitbar(0,waitText);
A = [];
b = [];
Aeq = [];
beq = [];
for i=1:maxIter
    if i>1
        psai0 = eqPoint_psai(:,i-1);
    else
        psai0 = initialPsai;
    end
    lb   = psai0 - 50*(pi/180);
    ub   = psai0 + 50*(pi/180);
    scale = force_field_symbolic(target_eqPoint_x(i),target_eqPoint_y(i),psai0,4);
    [eqPoint_psai(:,i),fval(i)] = fmincon(@(psai) force_field_symbolic(target_eqPoint_x(i),target_eqPoint_y(i),psai,4)/scale ,psai0,A,b,Aeq,beq,lb,ub);
    [all_eqPoints, fvals] = findEqPoints(all_eqPoints,x(2)-x(1),eqPoint_psai(:,i));
    all_eqPoints_x(:,i) = all_eqPoints(1,:)';
    all_eqPoints_y(:,i) = all_eqPoints(2,:)';
    waitbar(i/maxIter,waitHandle, waitText);
end
close(waitHandle)


% myPlotter(x,y,Frho,Faxial,MagPos,initialPsai, 0, eqPoints);
% myPlotter(x,y,Frho,Faxial,MagPos,newPsai, 0, eqPoints);


fid = fopen('psaiController.m', 'wt');
fprintf(fid, 'function [Psai,target_x,target_y,targetInd,real_x,real_y] = psaiController(t)\n');
fprintf(fid, 'dt = %f;\n', dt);
fprintf(fid, 'ind = floor(t/dt) + 1;\n');
fprintf(fid, 'ind = min(ind, %d);\n', maxIter);
fprintf(fid, 'target_x_sequence = [\n');
for i=1:size(target_eqPoint_x,1)
    fprintf(fid, '\t');
    fprintf(fid, '%f ', target_eqPoint_x(i,:));
    fprintf(fid, '\n');
end
fprintf(fid, '\n];\n');
fprintf(fid, 'target_y_sequence = [\n');
for i=1:size(target_eqPoint_y,1)
    fprintf(fid, '\t');
    fprintf(fid, '%f ', target_eqPoint_y(i,:));
    fprintf(fid, '\n');
end
fprintf(fid, '\n];\n');
fprintf(fid, 'Psai_sequence = [\n');
for i=1:size(eqPoint_psai,1)
    fprintf(fid, '\t');
    fprintf(fid, '%f ', eqPoint_psai(i,:));
    fprintf(fid, '\n');
end
fprintf(fid, '\n];\n');
fprintf(fid, 'real_x_sequence = [\n');
for i=1:size(all_eqPoints_x,1)
    fprintf(fid, '\t');
    fprintf(fid, '%f ', all_eqPoints_x(i,:));
    fprintf(fid, '\n');
end
fprintf(fid, '\n];\n');
fprintf(fid, 'real_y_sequence = [\n');
for i=1:size(all_eqPoints_y,1)
    fprintf(fid, '\t');
    fprintf(fid, '%f ', all_eqPoints_y(i,:));
    fprintf(fid, '\n');
end
fprintf(fid, '\n];\n');
fprintf(fid, 'Psai = Psai_sequence(:,ind);\n');
fprintf(fid, 'if nargout > 1\n');
fprintf(fid, '\ttargetInd = ind;\n');
fprintf(fid, '\ttarget_x = target_x_sequence;\n');
fprintf(fid, '\ttarget_y = target_y_sequence;\n');
fprintf(fid, '\treal_x = real_x_sequence(:,ind);\n');
fprintf(fid, '\treal_y = real_y_sequence(:,ind);\n');
fprintf(fid, 'end\n');
fprintf(fid, 'end');
fclose(fid);




% for i=1:Npoints
%     for j=1:Npoints
%         [a1,fval1(i,j)] = fminsearch(@(psai)abs(force_field_symbolic(x(i),y(j),psai,1)), [0,0,0,0]');
%         [a2,fval2(i,j)] = fminsearch(@(psai)abs(force_field_symbolic(x(i),y(j),psai,2)), [0,0,0,0]');
%         [a3,fval3(i,j)] = fminsearch(@(psai)abs(force_field_symbolic(x(i),y(j),psai,4)), [0,0,0,0]');
%     end
% end



function [eqPoints, fvals] = findEqPoints(initialEqPoints,resolution,Psai)
    eqPoints    = initialEqPoints;
    %
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    for i=1:size(initialEqPoints,2)
        lb = initialEqPoints(:,i)-5*resolution;
        ub = initialEqPoints(:,i)+5*resolution;
        ans0 = initialEqPoints(:,i);
        scale = force_field_symbolic(ans0(1),ans0(2),Psai,4);
        [ans1,fval(:,i)] = fmincon(@(x_) force_field_symbolic(x_(1),x_(2),Psai,4)/scale ,ans0,A,b,Aeq,beq,lb,ub);
        if fval(:,i) < 0.05
            eqPoints(:,i) = ans1;
        else
            eqPoints(:,i) = [-10,-10];
        end
    end
    if nargout > 1
        fvals = fval;
    end
end
%
%
function myPlotter(x,y,Frho,Faxial,MagPos,psai, showLocalMinLine, localMin)
    figHandle = figure;
    hold on
    if showLocalMinLine
        F = sqrt(Frho.^2+Faxial.^2);
        for i=1:length(x)
            [TF,P] = islocalmin(F(:,i),'FlatSelection','center');
            [row,~] = find(TF);
            col = ones(length(row))*i;
            plot(x(col), y(row), 'r*')
        end
        for i=1:length(x)
            [TF,P] = islocalmin(F(i,:),'FlatSelection','center');
            [~,col] = find(TF);
            row = ones(length(col))*i;
            plot(x(col), y(row), 'b*')
        end
    end
    if ~isempty(localMin)
        plot(localMin(1,:), localMin(2,:), 'b.', 'MarkerSize', 20);
    end
    plot_field = streamslice(x,y,Frho,Faxial,'method','cubic');
    set(plot_field,'Color','black','LineWidth',1.2);
    hold on 
    ax1 = plot(MagPos(:,1), MagPos(:,2), 'r.', 'MarkerSize', 15);
    title('\textbf{Force Field}','interpreter','latex', 'fontsize',14)
    xlabel('x [m]','interpreter','latex')
    ylabel('y [m]','interpreter','latex')
    xlim([min(x) max(x)]);
    ylim([min(y) max(y)]);
    axis square
    hold off
%     b = brush(figHandle);
%     b.Color = 'g';
%     b.Enable = 'on';

%     [TF,P] = islocalmin(sqrt(Frho.^2+Faxial.^2),'FlatSelection','center');
%     [row,col] = find(TF);
%     plot(x([col,row]), y([col,row]), 'r*')
    if ~isempty(psai)
        for i=1:length(psai)
%             titleText{i} = sprintf('${\\psi}_%d=%.0f^{\\circ}$  ', i, psai(i)*(180/pi));
            text{i} = sprintf('{\\psi}_%d=%.0f^{\\circ}\n', i, psai(i)*(180/pi));
        end
%         title(sprintf('%s', titleText{:}) ,'interpreter','latex', 'fontsize',14)
        annotation( 'textbox', 'String', sprintf('%s', text{:}), 'Color', 'k', ...
            'FontSize', 14, 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.75,0.9,0.5,0] )
%         set( gca, 'Position', [0.1, 0.1, 0.6, 0.8] )
    end
end