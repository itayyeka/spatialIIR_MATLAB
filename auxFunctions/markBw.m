function [bwVal] = markBw(azVec, bp, thetaS)
    bwAux_diff = abs(bp - db(1/2));
    [~,bwId] = min(bwAux_diff);
    bwVal = abs(azVec(bwId) - thetaS);
    bwVec = thetaS + sort([-1 1] * bwVal);
    for bw = bwVec
        xline(bw,'--r');
    end
    yline(db(1/2),'--r');
end