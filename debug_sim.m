function [] = debug_sim()
    clear all
    close all
    clc
    
    simCfg                  = [];
    simCfg.nSensors         = 3;
    simCfg.nIterations      = 100;
    simCfg.inputFreq        = 1e9;
    simCfg.nTheta           = 300;
    simCfg.r                = 0.7;
    simCfg.thetaS           = pi/2;    
    simCfg.lambdaD_ratio    = 1/2;
    
    %% run sim
    simOut                  = spatialIIR_singleFreq(simCfg);
    simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);
    theoryBp_dbAbs_norm     = f_theoryBp_dbAbs_norm(simOut.cfg.targetDuVec,simOut.cfg);
    figure;
    hold on;
    plot(simOut.cfg.targetAngleVec/pi,theoryBp_dbAbs_norm(:),'*b-');
    plot(simOut.cfg.targetAngleVec/pi,simBp_dbAbs_norm(:),'og:');
    
end