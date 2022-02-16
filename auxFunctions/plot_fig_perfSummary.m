function [] = plot_fig_perfSummary()
    simCfg              = [];
    simCfg.nSensors     = 20;
    simCfg.nIterations  = 100;
    simCfg.inputFreq    = 10e9;
    simCfg.nTheta       = 600;
    simCfg.thetaS               = pi/2;
    simCfg.targetRange_samples  = 20*simCfg.nSensors;
    
    simCfg_base     = simCfg;
    f_apertureImprove = @(r) 1.4/((1-r)*(1.4-0.4*r));
    f_expBw = @(N,r) 1.4/(f_apertureImprove(r)*N);
    
    if true
        figure;
        %% r=0
        simCfg                  = simCfg_base;
        simCfg.r                = 0;
        simOut                  = spatialIIR_singleFreq(simCfg);
        targetAngleVec          = simOut.cfg.targetAngleVec;
        targetDuVec             = simOut.cfg.targetDuVec;
        simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);
        subplot(2,2,1);
        plot(targetAngleVec/pi,simBp_dbAbs_norm(:));
        ylim([-50,0])
        bw00 = markBw(targetAngleVec/pi, simBp_dbAbs_norm(:),simOut.cfg.thetaS/pi)
        expBw00 = f_expBw(simCfg.nSensors, 0)
        err00 = abs(expBw00 - bw00)
        err00_prec = 100*err00/bw00
        hold on;
        %% r=0.3
        simCfg                  = simCfg_base;
        simCfg.r                = 0.3;
        simOut                  = spatialIIR_singleFreq(simCfg);
        targetAngleVec          = simOut.cfg.targetAngleVec;
        targetDuVec             = simOut.cfg.targetDuVec;
        simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);
        subplot(2,2,2);
        plot(targetAngleVec/pi,simBp_dbAbs_norm(:));
        ylim([-50,0])
        bw03 = markBw(targetAngleVec/pi, simBp_dbAbs_norm(:),simOut.cfg.thetaS/pi)
        expBw03 = f_expBw(simCfg.nSensors, 0.3)
        err03 = abs(expBw03 - bw03)
        err03_prec = 100*err03/bw03
        hold on;
        %% r=0.6
        simCfg                  = simCfg_base;
        simCfg.r                = 0.6;
        simOut                  = spatialIIR_singleFreq(simCfg);
        targetAngleVec          = simOut.cfg.targetAngleVec;
        targetDuVec             = simOut.cfg.targetDuVec;
        simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);
        subplot(2,2,3);
        plot(targetAngleVec/pi,simBp_dbAbs_norm(:));
        ylim([-50,0])
        bw06 = markBw(targetAngleVec/pi, simBp_dbAbs_norm(:),simOut.cfg.thetaS/pi)
        expBw06 = f_expBw(simCfg.nSensors, 0.6)
        err06 = abs(expBw06 - bw06)
        err06_prec = 100*err06/bw06
        hold on;
        %% r=0.8
        simCfg                  = simCfg_base;
        simCfg.r                = 0.8;
        simOut                  = spatialIIR_singleFreq(simCfg);
        targetAngleVec          = simOut.cfg.targetAngleVec;
        targetDuVec             = simOut.cfg.targetDuVec;
        simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);
        subplot(2,2,4);
        plot(targetAngleVec/pi,simBp_dbAbs_norm(:));
        ylim([-60,0])
        bw08 = markBw(targetAngleVec/pi, simBp_dbAbs_norm(:),simOut.cfg.thetaS/pi)
        expBw08 = f_expBw(simCfg.nSensors, 0.8)
        err08 = abs(expBw08 - bw08)
        err08_prec = 100*err08/bw08
        hold on;
    end    
    fixfig();
end