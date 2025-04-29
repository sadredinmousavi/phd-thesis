function [eqPoints] = convertEqPointRelativeTimeToAbs(eqPoints)


for i=2:length(eqPoints)
    eqPoints(1,i) = eqPoints(1,i-1) + eqPoints(1,i);
end


end

