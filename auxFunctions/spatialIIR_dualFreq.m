function [] = spatialIIR_dualFreq(simCfg)
try
    simCfg;
catch
    close all;
    simCfg                          = [];
    simCfg.nSensors                 = 5;
    simCfg.nIterations              = 100;
    simCfg.inputFreq                = 10e9;
    simCfg.dF_ratio                 = 1e-2;
    simCfg.r                        = 0.4;
    simCfg.rangeErrorToLambdaRatio  = 1000;
    simCfg.thetaS                   = (3/2)*pi/2;
    simCfg.snr                      = inf;
end
simCfg_base     = simCfg;
N               = simCfg.nSensors;

if true
    
    %% f1
    if true
        simCfg                      = simCfg_base;
        simCfg.compensationFreq     = 0*simCfg_base.inputFreq*simCfg.dF_ratio;
        simOut_f1                   = spatialIIR_singleFreq(simCfg);
        stft1                       = simOut_f1.stftMat(end,:);
    end
    %% f2
    if true
        simCfg                          = simCfg_base;
        simCfg.overrideFeedbackCoeffs   = [zeros(1,N-1) 1];
        simCfg.inputFreq                = simCfg_base.inputFreq*(1+simCfg.dF_ratio);
        simOut_f2                       = spatialIIR_singleFreq(simCfg);
        stft2                           = simOut_f2.stftMat(end,:);
        hTwoFreq_err0                   = 1./(1./stft1 - 1./stft2);
        hTwoFreq_err0_dbAbs2            = db(abs(hTwoFreq_err0));
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
        plot(targetAngleVec,hTwoFreq_err0_dbAbs2_norm(:),'*r-','MarkerIndices',3:4:nTheta);
    end
end
end