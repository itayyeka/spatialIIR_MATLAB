function [simOut] = spatialIIR_singleFreq(simCfg)
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
        r   = 0.95;
    end
    simCfg.r = r;
    try
        sigFreq     = simCfg.sigFreq;
    catch
        sigFreq     = 10e9;%Hz
    end
    simCfg.sigFreq = sigFreq;
    try
        compensationFreq     = simCfg.compensationFreq;
    catch
        compensationFreq     = sigFreq;
    end
    simCfg.compensationFreq = compensationFreq;
    try
        rangeErrorToLambdaRatio     = simCfg.rangeErrorToLambdaRatio;
    catch
        rangeErrorToLambdaRatio     = 0;
    end
    simCfg.rangeErrorToLambdaRatio = rangeErrorToLambdaRatio;    
    try
        snr = simCfg.snr;
    catch
        snr = inf;
    end
    try
        noiseType = simCfg.noiseType;
    catch
        noiseType = "output";%/"thermal"
    end
    simCfg.noiseType = noiseType;
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

    simOut.cfg = simCfg;

    %% auxiliary
    c               = propagationVelocity;
    lambda          = c/sigFreq;
    rangeError      = lambda*rangeErrorToLambdaRatio;
    try
        D   = simCfg.D;
    catch
        try
            D = lambda * simCfg.lambdaD_ratio;
        catch
            D   = lambda/2;
        end
    end
    simCfg.D        = D;
    f_exp           = @(f,t) exp(1i*2*pi*f*t);
    f_sig           = @(f,t) f_exp(f,t).*heaviside(t);
    randNorm        = (10^(-snr/10))/sqrt(2);
    f_noise_thermal = @(t) randNorm*exp(1i*2*pi*rand(size(t(:),1),nSensors));
    f_noise_output  = @(t) randNorm*exp(1i*2*pi*rand(size(t)));
    N               = nSensors;
    f_dTOA          = @(theta) reshape((0:(N-1))*D*cos(theta)/c,[],1);
    f_steering      = @(theta) fliplr(f_exp(sigFreq,f_dTOA(theta)));
    try
        simCfg.overrideAlpha;
        f_alphaVec  = @(thetaS,range) reshape(simCfg.overrideAlpha,[],1);
    catch
        f_alphaVec  = @(thetaS,range) ...
            reshape( ...
            (1/N) ...
            * ...
            f_steering(thetaS) ...
            * ...
            exp(-1i*2*pi*compensationFreq*2*(range+rangeError)/c) ...
            ,[],1);
    end
    try
        simCfg.overrideBeta;
        f_betaVec   = @(thetaS,range) reshape(simCfg.overrideBeta,[],1);
    catch
        f_betaVec   = @(thetaS,range) f_alphaVec(thetaS,range);
    end

    fSample                 = 5*sigFreq;
    tSample                 = 1/fSample;
    tPd                     = targetRange_samples*tSample;
    targetRange             = c*tPd;
    nSamplesIter            = floor(2*tPd/tSample);
    epochIterRatio          = 2;
    nSamplesEpoch           = floor(nSamplesIter/epochIterRatio);
    simDuration_epochs      = 1 + ceil(simDuration_iterations * epochIterRatio);

    try
        targetAngleVec  = simCfg.targetAngleVec;
    catch
        targetDuVec     = f_gen_duVec(simCfg);
        targetAngleVec  = f_convert_du2angle(targetDuVec,simCfg);
    end
    targetDuVec = f_convert_angle2du(targetAngleVec,simCfg);

    simDuration_samples     = simDuration_epochs*nSamplesEpoch;
    hMat                    = zeros(simDuration_samples,nTheta);
    stftVec                 = zeros(1,nTheta);

    alphaVec        = f_alphaVec(thetaS,targetRange);
    alphaVecT       = transpose(alphaVec);
    alphaVecH       = conj(alphaVecT);
    betaVec         = f_betaVec(thetaS,targetRange);
    betaVecT        = transpose(betaVec);
    betaVecH        = conj(betaVecT);
    sampleIdVec     = 1 : nSamplesEpoch;
    targetAngleId   = 0;
    for targetAngle = targetAngleVec
        targetAngleId   = targetAngleId + 1;
        dTOAVec         = f_dTOA(targetAngle);
        steerVec        = f_steering(targetAngle);
        feedbackSig     = zeros(simDuration_samples,1);
        arrayOutput     = zeros(simDuration_samples,1);
        for epochId = 1 : simDuration_epochs
            epochSampleIdVec                 = (epochId-1)*nSamplesEpoch + sampleIdVec;
            epochTVec                        = (epochSampleIdVec-1)*tSample;
            feedbackGenerationTime          = epochTVec - 2*tPd;
            feedbackGenerationTime_samples  = round(feedbackGenerationTime/tSample);
            historySampleIdVec              = ...
                round(...
                min(feedbackGenerationTime_samples)-historyBufferSize*nSamplesEpoch ...
                : ...
                max(feedbackGenerationTime_samples)+historyBufferSize*nSamplesEpoch ...
                );
            historyTVec                     = (historySampleIdVec-1)*tSample;
            dTOAMat                         = repmat(reshape(dTOAVec,1,[]),nSamplesEpoch,1);
            feedbackGenerationTimeMat       = repmat(feedbackGenerationTime(:),1,N)+dTOAMat;
            curEpochInput_sig               = f_sig(sigFreq,feedbackGenerationTimeMat);
            feedbackSig_valid               = feedbackSig(max(1,historySampleIdVec));
            curEpochInput_feedback          = interp1(historyTVec(:),feedbackSig_valid(:),feedbackGenerationTimeMat,'spline');

            if false
                figure;
                plot(historyTVec,imag(feedbackSig_valid))
                hold on;
                for antId = 1 : size(feedbackGenerationTimeMat,2)
                    plot(feedbackGenerationTimeMat(:,antId),imag(curEpochInput_feedback(:,antId)));
                end
                close all;
            end

            epochArrayInput                  = curEpochInput_sig + r*curEpochInput_feedback;
            if strcmpi(simCfg.noiseType,"thermal")
                epochArrayInput = epochArrayInput + f_noise_thermal(epochTVec(:));
            end
            epochArrayFeedback               = reshape(epochArrayInput*alphaVecH(:),[],1);
            feedbackSig(epochSampleIdVec)    = epochArrayFeedback(:);

            if false
                %% DEBUG
                figure;
                plot(abs(feedbackSig));
                close all;
            end

            epochArrayOutput                 = reshape(epochArrayInput*betaVecH(:),[],1);
            if strcmpi(simCfg.noiseType,"output")
                epochArrayOutput = epochArrayOutput + f_noise_output(epochTVec(:));
            end
            arrayOutput(epochSampleIdVec)    = epochArrayOutput(:);
        end

        hMat(:,targetAngleId)   = arrayOutput(:);
        tmpStft = stft(arrayOutput(:),fSample);
        [~,binId] = max(abs(tmpStft(:,end)));
        binVal = tmpStft(binId,end);
        stftVec(targetAngleId) = binVal;

    end

    simOut.cfg.nSamplesIter     = nSamplesIter;
    simOut.cfg.D                = D;
    simOut.cfg.lambda           = lambda;
    simOut.cfg.targetAngleVec   = targetAngleVec;
    simOut.cfg.targetDuVec      = targetDuVec;
    simOut.hMat                 = hMat;
    simOut.stftVec              = stftVec;
end