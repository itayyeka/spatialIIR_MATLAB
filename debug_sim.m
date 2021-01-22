function [] = debug_sim()
    clear all
    close all
    clc
    
    simCfg              = [];
    simCfg.nSensors     = 3;
    simCfg.nIterations  = 100;
    simCfg.inputFreq    = 1e9;
    simCfg.nTheta       = 300;
    simCfg.r            = 0.7;
    simCfg.thetaS       = pi/3;
    
    defaultSimOut   = spatialIIR_singleFreq(simCfg);
    targetAngleVec  = defaultSimOut.targetAngleVec;
    N               = defaultSimOut.cfg.nSensors;
    omega           = 2*pi*simCfg.inputFreq;
    duCoef          = omega*defaultSimOut.cfg.D/defaultSimOut.cfg.propagationVelocity;
    duVec           = duCoef*(cos(targetAngleVec)-cos(defaultSimOut.cfg.thetaS));
    duVec_angle     = acos(cos(defaultSimOut.cfg.thetaS) + defaultSimOut.cfg.propagationVelocity*duVec/(omega*defaultSimOut.cfg.D));
    %% run sim
    theoryBp_dbAbs_norm     = f_theoryBp_dbAbs_norm(duVec,simCfg.r,N);
    simOut                  = spatialIIR_singleFreq(simCfg);
    simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);
    figure;
    hold on;
    plot(duVec_angle,theoryBp_dbAbs_norm(:),'*b-');
    plot(duVec_angle,simBp_dbAbs_norm(:),'og:');
    
end