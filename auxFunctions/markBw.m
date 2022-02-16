function [bwVal] = markBw(azVec_norm, bp, thetaS_norm)    
    az_interp   = linspace(min(azVec_norm),max(azVec_norm),10000);
    bp_interp   = interp1(azVec_norm,bp,az_interp);
    %% part1
    azVec_norm1 = az_interp(az_interp<=thetaS_norm);
    bp1         = bp_interp(az_interp<=thetaS_norm);
    bwAux_diff  = abs(bp1 - db(1/4));
    [~,bwId]    = min(bwAux_diff);
    bwVal1      = thetaS_norm - abs(azVec_norm1(bwId) - thetaS_norm);
    %% part2
    azVec_norm2 = az_interp(az_interp>thetaS_norm);
    bp1         = bp_interp(az_interp>thetaS_norm);
    bwAux_diff  = abs(bp1 - db(1/4));
    [~,bwId]    = min(bwAux_diff);
    bwVal2      = thetaS_norm + abs(azVec_norm2(bwId) - thetaS_norm);
    %% combine
    bwVec = [bwVal1 , bwVal2];
    bwVal = pi*(bwVal2 - bwVal1)/2;
    for bw = bwVec
        xline(bw,'--r');
    end
    yline(db(1/4),'--r');
end