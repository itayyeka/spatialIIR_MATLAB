function [duVec] = f_convert_angle2du(angleVec,simCfg)
    omega   = 2*pi*simCfg.sigFreq;
    duCoef  = (omega*simCfg.D)/simCfg.propagationVelocity;
    duVec   = duCoef * (cos(angleVec) - cos(simCfg.thetaS));
end

