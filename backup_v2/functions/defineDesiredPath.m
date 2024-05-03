function [desiredPath] = defineDesiredPath(point1, point2, startTime, endTime, dt)
    time = startTime:dt:endTime;
    for i=1:length(time)
        desiredPath(1,i) = time(i);
        desiredPath(2,i) = point1(1) + (point2(1)-point1(1)) * (time(i)-startTime)/(endTime-startTime);
        desiredPath(3,i) = point1(2) + (point2(2)-point1(2)) * (time(i)-startTime)/(endTime-startTime);
        %x_desiredPath(i) = x_desiredPath(1) + maxDisp*(time);
        %y_desiredPath(i) = -maxDisp*(time)^4 + y_desiredPath(1);
    end
end

