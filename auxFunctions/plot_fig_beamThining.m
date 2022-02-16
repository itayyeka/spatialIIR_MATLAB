function [] = plot_fig_beamThining()
    simCfg              = [];
    simCfg.nSensors     = 3;
    simCfg.nIterations  = 1;
    simCfg.sigFreq      = 10e9;
    simCfg.nTheta       = 300;
    simCfg.r            = 0.9;
    
    simCfg_base     = simCfg;
    defaultSimOut   = spatialIIR_singleFreq(simCfg);
    targetAngleVec  = defaultSimOut.cfg.targetAngleVec;
    nIterVec        = [1, 2, 5, 10, 50];
    figure;
    for nIter = nIterVec
        %% sim
        simCfg                  = simCfg_base;
        simCfg.nIterations      = nIter+1;
        simOut                  = spatialIIR_singleFreq(simCfg);
        simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);
        hold on;
        plot(targetAngleVec/pi,simBp_dbAbs_norm(:));
        %% cfg
%         ylim([-50 0]);
%         fixfig;
    end
end