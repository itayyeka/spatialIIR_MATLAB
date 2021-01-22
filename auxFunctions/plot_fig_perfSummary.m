function [] = plot_fig_perfSummary()
    simCfg              = [];
    simCfg.nSensors     = 3;
    simCfg.nIterations  = 100;
    simCfg.inputFreq    = 10e9;
    simCfg.nTheta       = 300;
    
    simCfg_base     = simCfg;
    defaultSimOut   = spatialIIR_singleFreq(simCfg);
    targetAngleVec  = defaultSimOut.targetAngleVec;
    if true
        figure;
        %% r=0
        simCfg                  = simCfg_base;
        simCfg.r                = 0;
        simOut                  = spatialIIR_singleFreq(simCfg);
        simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);
        subplot(2,2,1);
        plot(targetAngleVec/pi,simBp_dbAbs_norm(:));
        markBw(targetAngleVec/pi, simBp_dbAbs_norm(:),simOut.cfg.thetaS/pi)
        hold on;
        %% r=0.3
        simCfg                  = simCfg_base;
        simCfg.r                = 0.3;
        simOut                  = spatialIIR_singleFreq(simCfg);
        simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);
        subplot(2,2,2);
        plot(targetAngleVec/pi,simBp_dbAbs_norm(:));
        markBw(targetAngleVec/pi, simBp_dbAbs_norm(:),simOut.cfg.thetaS/pi)
        hold on;
        %% r=0.6
        simCfg                  = simCfg_base;
        simCfg.r                = 0.6;
        simOut                  = spatialIIR_singleFreq(simCfg);
        simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);
        subplot(2,2,3);
        plot(targetAngleVec/pi,simBp_dbAbs_norm(:));
        markBw(targetAngleVec/pi, simBp_dbAbs_norm(:),simOut.cfg.thetaS/pi)
        hold on;
        %% r=0.8
        simCfg                  = simCfg_base;
        simCfg.r                = 0.8;
        simOut                  = spatialIIR_singleFreq(simCfg);
        simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);
        subplot(2,2,4);
        plot(targetAngleVec/pi,simBp_dbAbs_norm(:));
        markBw(targetAngleVec/pi, simBp_dbAbs_norm(:),simOut.cfg.thetaS/pi)
        hold on;
    end
    fixfig();
end