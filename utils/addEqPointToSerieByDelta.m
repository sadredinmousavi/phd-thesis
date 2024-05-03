function [eqPoints] = addEqPointToSerieByDelta(xDelta, yDelta, timeDelta, eqPoints)
%ADDEQPOINTTOSERIELINEAR Summary of this function goes here

if nargin > 3
    lastData = eqPoints(:,end) + [timeDelta; xDelta; yDelta];
    eqPoints = [eqPoints lastData];
else
    eqPoints = [timeDelta; xDelta; yDelta];
end

end

