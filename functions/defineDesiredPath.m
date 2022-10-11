function [x_desiredPath,y_desiredPath] = defineDesiredPath(point,inputs)

    pointIndexInsideEqPoints = inputs.indexInsideEqPoints;
    dt = inputs.dt;
    startTime = inputs.startTime;
    endTime = inputs.endTime;
    totalTime = endTime - startTime;
    maxDisp = inputs.maxDisp;
    maxIter = endTime/dt + 1;
    x_desiredPath(1) = point(1);
    y_desiredPath(1) = point(2);
    for i=2:maxIter
        if dt*i >= startTime
            time = (dt*(i-1) - startTime) / totalTime; % from 0 to 1
            x_desiredPath(i) = x_desiredPath(1) + maxDisp*(time);
            y_desiredPath(i) = -maxDisp*(time)^4 + y_desiredPath(1);
        else
            x_desiredPath(i) = x_desiredPath(1);
            x_desiredPath(i) = x_desiredPath(1);
        end
    end
end

