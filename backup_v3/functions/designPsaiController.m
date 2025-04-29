function [outputs] = designPsaiController(configs, handles)

    time_tmp = configs.eqPoints{1}(1,:);
    for i=2:length(configs.eqPoints)
        time_tmp = [time_tmp configs.eqPoints{i}(1,:)];
    end
    time  = unique(time_tmp);
    eqP = {};
    lambda = {};
    %
    waitText  = 'openloop Psai calculations - Please wait...';
    waitHandle = waitbar(0,waitText);
    maxIter = length(time);
    
    for i=1:maxIter
        hasLambda = 0;
        for cnt = 1:length(configs.eqPoints)
            temp = configs.eqPoints{cnt}(1,:) - time(i);
            for k = 1:length(temp)
                if temp(k) > 0
                    eqP{cnt}(2:3,i) = configs.eqPoints{cnt}(2:3,k-1);
                    if size(configs.eqPoints{cnt}, 1) > 3
                        lambda{cnt} = configs.eqPoints{cnt}(4:end,k-1);
                        hasLambda = 1;
                    end
                    break;
                elseif k == length(temp)
                    eqP{cnt}(2:3,i) = configs.eqPoints{cnt}(2:3,end);
                    if size(configs.eqPoints{cnt}, 1) > 3
                        lambda{cnt} = configs.eqPoints{cnt}(4:end,end);
                        hasLambda = 1;
                    end
                end
            end
            eqP{cnt}(1,i) = time(i);
        end
        %
        %
        if configs.usePreparedPsai
            Psai_t(:,i) = configs.Psai(:,i);
        else
            for j=1:length(eqP)
                eqP_tmp{j}(:,1) = eqP{j}(2:3,i);
            end
            if isequal(handles.calcPsaiFromEqFunc, @calculatePsai_minimization)
                if eqP_tmp{1}(1) == eqP_tmp{2}(1) && eqP_tmp{1}(2) == eqP_tmp{2}(2)
%                     handles.calcPsaiFromEqFunc = @calculatePsai_minimization_2;
                end
            end
            if hasLambda
                [rankM, error, hasAns, isStable, Psai_t(:,i)] = handles.calcPsaiFromEqFunc(eqP_tmp,  configs.epms.MagPos, lambda);
            else
                [rankM, error, hasAns, isStable, Psai_t(:,i)] = handles.calcPsaiFromEqFunc(eqP_tmp,  configs.epms.MagPos);
            end
        end
        
        %
        if handles.findEqFromMinimization
%             [x_eqPoints(:,i),y_eqPoints(:,i)] = findEqPoints_Minimization(inputs.x_space,inputs.y_space,inputs.Psai_0);
            try
%                 [x_eqPoints(:,i),y_eqPoints(:,i)] = findEqPoints_Minimization(configs.x_space,configs.y_space,Psai_t(:,i));
                [x_eqPoints(:,i),y_eqPoints(:,i)] = findEqPoints_Dynamics(configs.vars.x_space,configs.vars.y_space,Psai_t(:,i));
            catch
%                 [x_eq_,y_eq_] = findEqPoints_Minimization(inputs.x_space,inputs.y_space,Psai_t(:,i));
%                 a = length(x_eqPoints(:,1)) - length(x_eq_);
%                 x_eqPoints(:,i) = [x_eq_ zeros(1,a)]';
%                 y_eqPoints(:,i) = [y_eq_ zeros(1,a)]';
                x_eqPoints(:,i) = x_eqPoints(:,1);
                y_eqPoints(:,i) = y_eqPoints(:,1);
            end
        else
            for cnt = 1:length(configs.eqPoints)
                x_eqPoints(cnt,i) = configs.eqPoints{cnt}(2,i);
                y_eqPoints(cnt,i) = configs.eqPoints{cnt}(3,i);
            end
        end
        eqPoints{i}.time = time(i);
        eqPoints{i}.x = x_eqPoints(:,i);
        eqPoints{i}.y = y_eqPoints(:,i);
        waitbar(i/maxIter,waitHandle, waitText);
    end
    close(waitHandle)

    outputs.time = time;
    outputs.eqP = eqP;
    outputs.Psai_t = Psai_t;
    outputs.x_eqPoints = x_eqPoints;
    outputs.y_eqPoints = y_eqPoints;





end

