function [] = spatialIIR_dualFreq(simCfg)
try
    simCfg;
catch
    close all;
    simCfg                          = [];
    simCfg.nSensors                 = 5;
    simCfg.nIterations              = 100;
    simCfg.inputFreq                = 10e9;
    simCfg.dF_ratio                 = 1e-1;
    simCfg.targetRange_samples      = 32;
    simCfg.r                        = 0.6;
    simCfg.rangeErrorToLambdaRatio  = 0.3;
    simCfg.thetaS                   = (1.5)*pi/2;
    simCfg.snr                      = inf;
    simCfg.kappa                    = simCfg.r;
end
simCfg_base     = simCfg;
N               = simCfg.nSensors;

if true
    %% singleFreq
    if true
        simCfg              = simCfg_base;
        simOut_singleReq    = spatialIIR_singleFreq(simCfg);
        bp_singleReq        = db(abs(simOut_singleReq .hMat(end,:)).^2);
        bp_singleReq_norm   = bp_singleReq-max(bp_singleReq(:));
    end
    %% f1
    if true
        simCfg                      = simCfg_base;
        simCfg.r                    = simCfg_base.r*simCfg_base.kappa;
        simCfg.compensationFreq     = 0*simCfg.dF_ratio*simCfg.inputFreq;
        simOut_f1                   = spatialIIR_singleFreq(simCfg);
        bp_f1                       = db(abs(simOut_f1.hMat(end,:)).^2);
        stft1                       = simCfg.r*simOut_f1.stftMat(end,:);
    end
    %% f2
    if true
        simCfg                          = simCfg_base;
        simCfg.overrideFeedbackCoeffs   = [zeros(1,N-1) 1];
        simCfg.inputFreq                = simCfg_base.inputFreq*(1+simCfg.dF_ratio);
        simOut_f2                       = spatialIIR_singleFreq(simCfg);
        bp_f2                           = db(abs(simOut_f2.hMat(end,:)).^2);
        stft2                           = simCfg.r*simOut_f2.stftMat(end,:);
        hTwoFreq_err0                   = 1./(1./stft1 - 1./stft2);
        hTwoFreq_err0_dbAbs2            = db(abs(hTwoFreq_err0).^2);
        hTwoFreq_err0_dbAbs2_norm       = hTwoFreq_err0_dbAbs2 - max(hTwoFreq_err0_dbAbs2(:));
    end
    %% theory
    targetAngleVec          = simOut_f1.targetAngleVec;
    duVec                   = pi*(cos(targetAngleVec)-cos(simCfg.thetaS));
    theoryBp_dbAbs_norm     = f_theoryBp_dbAbs_norm(duVec,simCfg_base.r,N);
    nTheta                  = simOut_f1.cfg.nTheta;
    %% plot
    if true
        figure;
        hold on;
        plot(targetAngleVec,theoryBp_dbAbs_norm(:),'*b-','MarkerIndices',1:4:nTheta);
        plot(targetAngleVec,bp_singleReq_norm,'*r-','MarkerIndices',2:4:nTheta);
        plot(targetAngleVec,hTwoFreq_err0_dbAbs2_norm(:),'*g-','MarkerIndices',3:4:nTheta);
        legend({'theory' 'singleFreq' 'dualFreq'});
    end
end
end