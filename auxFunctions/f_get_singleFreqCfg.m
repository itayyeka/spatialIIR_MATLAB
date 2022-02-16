function [simCfg] = f_get_singleFreqCfg(simCfg)
    %% cfg
    try
        simDuration_iterations  = simCfg.nIterations;
    catch
        simDuration_iterations  = 100;%sec
    end
    simCfg.simDuration_iterations = simDuration_iterations;
    try
        nSensors    = simCfg.nSensors;
    catch
        nSensors    = 3;
    end
    simCfg.nSensors = nSensors;
    try
        r   = simCfg.r;
    catch
        r   = 0.6;
    end
    simCfg.r = r;
    try
        sigFreq     = simCfg.sigFreq;
    catch
        sigFreq     = 10e9;%Hz
    end
    simCfg.sigFreq = sigFreq;    
    try
        rangeErrorToLambdaRatio     = simCfg.rangeErrorToLambdaRatio;
    catch
        rangeErrorToLambdaRatio     = 0;
    end
    simCfg.rangeErrorToLambdaRatio = rangeErrorToLambdaRatio;
    try
        stftDuration_iterPrecent    = simCfg.stftDuration_iterPrecent;
    catch
        stftDuration_iterPrecent    = 90;
    end
    simCfg.stftDuration_iterPrecent = stftDuration_iterPrecent;
    try
        snr = simCfg.snr;
    catch
        snr = inf;
    end
    simCfg.snr = snr;
    try
        thetaS  = simCfg.thetaS;
    catch
        thetaS  = pi/2;
    end
    simCfg.thetaS = thetaS;
    try
        targetRange_samples = simCfg.targetRange_samples;
    catch
        targetRange_samples = 150;
    end
    simCfg.targetRange_samples = targetRange_samples;
    try
        nTheta              = simCfg.nTheta;
    catch
        nTheta              = 100;
    end
    simCfg.nTheta = nTheta;
    try
        thetaMin            = simCfg.thetaMin;
        thetaMax            = simCfg.thetaMax;
        simCfg.thetaMin     = thetaMin;
        simCfg.thetaMax     = thetaMax;
    catch
    end
    try
        duMin           = simCfg.duMin;
        duMax           = simCfg.duMax;
        simCfg.duMin    = duMin;
        simCfg.duMax    = duMax;
    catch
    end
    try
        propagationVelocity    = simCfg.propagationVelocity;
    catch
        propagationVelocity    = 3e8;
    end
    simCfg.propagationVelocity = propagationVelocity;
    try
        historyBufferSize    = simCfg.historyBufferSize;
    catch
        historyBufferSize    = 0.7;
    end
    simCfg.historyBufferSize = historyBufferSize;
       
    %% auxiliary
    c               = propagationVelocity;
    lambda          = c/sigFreq;
    try
        D   = simCfg.D;
    catch
        try
            D = lambda * simCfg.lambdaD_ratio;
        catch
            D   = lambda/2;
        end
    end
    simCfg.D = D;
end