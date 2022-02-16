function [angleVec] = f_convert_du2angle(duVec,simCfg)
    omega       = 2*pi*simCfg.sigFreq;
    duCoef      = simCfg.propagationVelocity/(omega*simCfg.D);
    angleVec    = acos(duCoef*duVec + cos(simCfg.thetaS));
end

