function [] = plot_fig_bpPresent()
    close all;
    clear all;
    clc;
    
    simCfg              = [];
    simCfg.r            = 0;
    simCfg.nIterations  = 100;
    simCfg.inputFreq    = 1e9;
    simCfg.nTheta       = 200;    
    simCfg_base         = simCfg;    
    
    %% Compute Bp    
    NVec    = [3, 6, 12];
    
    nBp         = numel(NVec);
    angleMat    = zeros(simCfg.nTheta,nBp);
    duMat       = zeros(simCfg.nTheta,nBp);
    bpMat       = zeros(simCfg.nTheta,nBp);
    for bpId = 1:nBp
        simCfg              = simCfg_base;
        simCfg.nSensors     = NVec(bpId);
        simOut              = spatialIIR_singleFreq(simCfg);
        angleMat(:,bpId)    = simOut.cfg.targetAngleVec(:);
        duMat(:,bpId)       = f_convert_angle2du(angleMat(:,bpId),simOut.cfg);
        bpMat(:,bpId)       = calc_simBp_dbAbs_norm(simOut);
    end
    figure;
    plot(angleMat,bpMat);
    ylim([-70, 0])
    fixfig;
    figure;
    polardb(angleMat,bpMat,-60);    
end