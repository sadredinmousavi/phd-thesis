function [] = designPsaiController(inputs)

    eqPoint1 = inputs.eqPoint1;
    eqPoint2 = inputs.eqPoint2;
    time1 = eqPoint1(1,:);
    time2 = eqPoint2(1,:);
    time  = unique([time1 time2]);
    %
    waitText  = 'openloop Psai calculations - Please wait...';
    waitHandle = waitbar(0,waitText);
    maxIter = length(time);
    for i=1:maxIter
        temp = time1 - time(i);
        for k = 1:length(temp)
            if temp(k) > 0
                eqP1(:,i) = eqPoint1(:,k-1);
                break;
            elseif k == length(temp)
                eqP1(:,i) = eqPoint1(:,end);
            end
        end
        eqP1(1,i) = time(i);
        %
        temp = time2 - time(i);
        for k = 1:length(temp)
            if temp(k) > 0
                eqP2(:,i) = eqPoint2(:,k-1);
                break;
            elseif k == length(temp)
                eqP2(:,i) = eqPoint2(:,end);
            end
        end
        eqP2(1,i) = time(i);
        %
        if inputs.usePrepaidPsai
            Psai_t(:,i) = inputs.Psai(:,i);
        else
            [rankM, error, hasAns, isStable, Psai_t(:,i)] = inputs.calcPsaiFromEqFunc([eqP1(2,i), eqP1(3,i)], [eqP2(2,i), eqP2(3,i)],  inputs.MagPos);
        end
        
        %
        if inputs.findEqFromMinimization
%             [x_eqPoints(:,i),y_eqPoints(:,i)] = findEqPoints_Minimization(inputs.x_space,inputs.y_space,inputs.Psai_0);
            try 
                [x_eqPoints(:,i),y_eqPoints(:,i)] = findEqPoints_Minimization(inputs.x_space,inputs.y_space,Psai_t(:,i));
            catch
%                 [x_eq_,y_eq_] = findEqPoints_Minimization(inputs.x_space,inputs.y_space,Psai_t(:,i));
%                 a = length(x_eqPoints(:,1)) - length(x_eq_);
%                 x_eqPoints(:,i) = [x_eq_ zeros(1,a)]';
%                 y_eqPoints(:,i) = [y_eq_ zeros(1,a)]';
                x_eqPoints(:,i) = x_eqPoints(:,1);
                y_eqPoints(:,i) = y_eqPoints(:,1);
            end
        else
            x_eqPoints(:,i) = [eqP1(2,i); eqP2(2,i)];
            y_eqPoints(:,i) = [eqP1(3,i); eqP2(3,i)];
        end
        eqPoints{i}.time = time(i);
        eqPoints{i}.x = x_eqPoints(:,i);
        eqPoints{i}.y = y_eqPoints(:,i);
        waitbar(i/maxIter,waitHandle, waitText);
    end
%     for i=1:maxIter
%         if i < maxIter/3
%             [Psai_t(:,i)] = ones(1,size(inputs.MagPos,1))* 090 *pi/180;
%         else
%             [Psai_t(:,i)] = (pi/180) * [20; 20; 60]; % ones(1,size(inputs.MagPos,1))* 040 *pi/180;
%         end
%         x_eqPoints(:,i) = 0;
%         y_eqPoints(:,i) = 0;
%     end
    close(waitHandle)



    fid = fopen('psaiController.m', 'wt');
    fprintf(fid, 'function [Psai,eqPoint1,eqPoint2,targetInd,real_eqPoint_x,real_eqPoint_y] = psaiController(t)\n');
    fprintf(fid, 'time = [');
    fprintf(fid, '%f ', time);
    fprintf(fid, '];\n');
    fprintf(fid, 'temp = time - t;\n');
    fprintf(fid, 'for k = 1:length(temp)\n');
    fprintf(fid, '\tif temp(k) > 0\n');
    fprintf(fid, '\t\tind = k-1;\n');
    fprintf(fid, '\t\tbreak;\n');
    fprintf(fid, '\telseif k == length(temp)\n');
    fprintf(fid, '\t\tind = k;\n');
    fprintf(fid, '\tend\n');
    fprintf(fid, 'end\n');
    fprintf(fid, 'eqPoint1 = [\n');
    for i=2:3
        fprintf(fid, '\t');
        fprintf(fid, '%f ', eqP1(i,:));
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n];\n');
    fprintf(fid, 'eqPoint2 = [\n');
    for i=2:3
        fprintf(fid, '\t');
        fprintf(fid, '%f ', eqP2(i,:));
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
    fprintf(fid, '\teqPoint1 = eqPoint1(:,ind);\n');
    fprintf(fid, '\teqPoint2 = eqPoint2(:,ind);\n');
    fprintf(fid, '\treal_eqPoint_x = real_eqPoint_x_sequence(:,ind);\n');
    fprintf(fid, '\treal_eqPoint_y = real_eqPoint_y_sequence(:,ind);\n');
    fprintf(fid, 'end\n');
    fprintf(fid, 'end');
    fclose(fid);



end

