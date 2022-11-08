function [] = designPsaiController(inputs)

    time_tmp = inputs.eqPoints{1}(1,:);
    for i=2:length(inputs.eqPoints)
        time_tmp = [time_tmp inputs.eqPoints{i}(1,:)];
    end
    time  = unique(time_tmp); 
    %
    waitText  = 'openloop Psai calculations - Please wait...';
    waitHandle = waitbar(0,waitText);
    maxIter = length(time);
    
    for i=1:maxIter
        for cnt = 1:length(inputs.eqPoints)
            temp = inputs.eqPoints{cnt}(1,:) - time(i);
            for k = 1:length(temp)
                if temp(k) > 0
                    eqP{cnt}(:,i) = inputs.eqPoints{cnt}(:,k-1);
                    break;
                elseif k == length(temp)
                    eqP{cnt}(:,i) = inputs.eqPoints{cnt}(:,end);
                end
            end
            eqP{cnt}(1,i) = time(i);
        end
        %
        %
        if inputs.usePrepaidPsai
            Psai_t(:,i) = inputs.Psai(:,i);
        else
            [rankM, error, hasAns, isStable, Psai_t(:,i)] = inputs.calcPsaiFromEqFunc(eqP,  inputs.MagPos);
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
            for cnt = 1:length(inputs.eqPoints)
                x_eqPoints(cnt,i) = inputs.eqPoints{cnt}(2,i);
                y_eqPoints(cnt,i) = inputs.eqPoints{cnt}(3,i);
            end
        end
        eqPoints{i}.time = time(i);
        eqPoints{i}.x = x_eqPoints(:,i);
        eqPoints{i}.y = y_eqPoints(:,i);
        waitbar(i/maxIter,waitHandle, waitText);
    end
    close(waitHandle)



    fid = fopen('psaiController.m', 'wt');
    fprintf(fid, 'function [Psai,eqPoints,targetInd,real_eqPoint_x,real_eqPoint_y] = psaiController(t)\n');
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
    for i=1:length(eqP)
        fprintf(fid, 'eqPoints_sequence{%d} = [\n', i);
        for k=2:3
            fprintf(fid, '\t');
            fprintf(fid, '%f ', eqP{i}(k,:));
            fprintf(fid, '\n');
        end
        fprintf(fid, '\n];\n');
    end
    fprintf(fid, 'Psai_sequence = [\n');
    for i=1:size(Psai_t,1)
        fprintf(fid, '\t');
        fprintf(fid, '%f ', Psai_t(i,:));
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n];\n');
    fprintf(fid, 'real_eqPoints_x_sequence = [\n');
    for i=1:size(x_eqPoints,1)
        fprintf(fid, '\t');
        fprintf(fid, '%f ', x_eqPoints(i,:));
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n];\n');
    fprintf(fid, 'real_eqPoints_y_sequence = [\n');
    for i=1:size(y_eqPoints,1)
        fprintf(fid, '\t');
        fprintf(fid, '%f ', y_eqPoints(i,:));
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n];\n');
    fprintf(fid, 'Psai = Psai_sequence(:,ind);\n');
    fprintf(fid, 'if nargout > 1\n');
    fprintf(fid, '\ttargetInd = ind;\n');
%     fprintf(fid, '\teqPoints = eqPoints(:,ind);\n');
    fprintf(fid, '\teqPoints = eqPoints_sequence;\n');
    fprintf(fid, '\treal_eqPoint_x = real_eqPoints_x_sequence(:,ind);\n');
    fprintf(fid, '\treal_eqPoint_y = real_eqPoints_y_sequence(:,ind);\n');
    fprintf(fid, 'end\n');
    fprintf(fid, 'end');
    fclose(fid);



end

