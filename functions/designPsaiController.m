function [x_eqPoints_0] = designPsaiController(inputs)

    x_eqPoints_0 = inputs.eqPoints.x_eqPoints_0;
    y_eqPoints_0 = inputs.eqPoints.y_eqPoints_0;
    x_desiredPath1 = inputs.dPath.x_desiredPath1;
    y_desiredPath1 = inputs.dPath.y_desiredPath1;
    x_desiredPath2 = inputs.dPath.x_desiredPath2;
    y_desiredPath2 = inputs.dPath.y_desiredPath2;
    %
    indexInsideEqPoints = inputs.dPath.indexInsideEqPoints;
    dt = inputs.dPath.dt;
    endTime = inputs.dPath.endTime;
    maxDisp = inputs.dPath.maxDisp;
    maxIter = endTime/dt + 1;
    %
    %

    waitText  = 'openloop Psai calculations - Please wait...';
    waitHandle = waitbar(0,waitText);
%     for i=1:maxIter       
% %         Psai_t(:,i) = inputs.Psai_0;
%         [Psai_t(:,i)] = calculatePsai_3PM([x_desiredPath1(i), y_desiredPath1(i)], inputs.MagPos, 30*pi/180)';
% %         [Psai_t(:,i)] = calculatePsai_6PM([x_desiredPath1(i), y_desiredPath1(i)], [x_desiredPath2(i), y_desiredPath2(i)],  inputs.MagPos)';
%         [x_eqPoints(:,i),y_eqPoints(:,i)] = findEqPoints_Minimization(inputs.x_space,inputs.y_space,inputs.Psai_0);
% %         try 
% %             [x_eqPoints(:,i),y_eqPoints(:,i)] = findEqPoints_Minimization(inputs.x_space,inputs.y_space,Psai_t(:,i));
% %         catch
% %             [x_eq_,y_eq_] = findEqPoints_Minimization(inputs.x_space,inputs.y_space,Psai_t(:,i));
% %             a = length(x_eqPoints(:,1)) - length(x_eq_);
% %             x_eqPoints(:,i) = [x_eq_ zeros(1,a)]';
% %             y_eqPoints(:,i) = [y_eq_ zeros(1,a)]';
% %         end        
%         waitbar(i/maxIter,waitHandle, waitText);
%     end
    for i=1:maxIter
        if i < maxIter/3
            [Psai_t(:,i)] = ones(1,size(inputs.MagPos,1))* 090 *pi/180;
        else
            [Psai_t(:,i)] = (pi/180) * [20; 20; 60]; % ones(1,size(inputs.MagPos,1))* 040 *pi/180;
        end
        x_eqPoints(:,i) = 0;
        y_eqPoints(:,i) = 0;
    end
    close(waitHandle)



    fid = fopen('psaiController.m', 'wt');
    fprintf(fid, 'function [Psai,target_x,target_y,targetInd,real_eqPoint_x,real_eqPoint_y] = psaiController(t)\n');
    fprintf(fid, 'dt = %f;\n', dt);
    fprintf(fid, 'ind = floor(t/dt) + 1;\n');
    fprintf(fid, 'ind = min(ind, %d);\n', maxIter);
    fprintf(fid, 'target_x_sequence = [\n');
    for i=1:size(x_desiredPath1,1)
        fprintf(fid, '\t');
        fprintf(fid, '%f ', x_desiredPath1(i,:));
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n];\n');
    fprintf(fid, 'target_y_sequence = [\n');
    for i=1:size(y_desiredPath1,1)
        fprintf(fid, '\t');
        fprintf(fid, '%f ', y_desiredPath1(i,:));
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n];\n');
    fprintf(fid, 'Psai_sequence = [\n');
    for i=1:size(Psai_t,1)
        fprintf(fid, '\t');
        fprintf(fid, '%f ', Psai_t(i,:));
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n];\n');
    fprintf(fid, 'real_eqPoint_x_sequence = [\n');
    for i=1:size(x_eqPoints,1)
        fprintf(fid, '\t');
        fprintf(fid, '%f ', x_eqPoints(i,:));
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n];\n');
    fprintf(fid, 'real_eqPoint_y_sequence = [\n');
    for i=1:size(y_eqPoints,1)
        fprintf(fid, '\t');
        fprintf(fid, '%f ', y_eqPoints(i,:));
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n];\n');
    fprintf(fid, 'Psai = Psai_sequence(:,ind);\n');
    fprintf(fid, 'if nargout > 1\n');
    fprintf(fid, '\ttargetInd = ind;\n');
    fprintf(fid, '\ttarget_x = target_x_sequence;\n');
    fprintf(fid, '\ttarget_y = target_y_sequence;\n');
    fprintf(fid, '\treal_eqPoint_x = real_eqPoint_x_sequence(:,ind);\n');
    fprintf(fid, '\treal_eqPoint_y = real_eqPoint_y_sequence(:,ind);\n');
    fprintf(fid, 'end\n');
    fprintf(fid, 'end');
    fclose(fid);



end

