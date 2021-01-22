function [simOut] = spatialIIR_singleFreq(simCfg)
    %% cfg
    
    try
        simDuration_iterations  = simCfg.nIterations;
    catch
        simDuration_iterations  = 100;%sec
    end
    try
        nSensors    = simCfg.nSensors;
    catch
        nSensors    = 3;
    end
    try
        r   = simCfg.r;
    catch
        r   = 0.95;
    end
    try
        sigFreq     = simCfg.inputFreq;
    catch
        sigFreq     = 10e9;%Hz
    end
    try
        compensationFreq    = simCfg.compensationFreq;
    catch
        compensationFreq    = sigFreq;
    end
    try
        rangeErrorToLambdaRatio     = simCfg.rangeErrorToLambdaRatio;
    catch
        rangeErrorToLambdaRatio     = 0;
    end
    try
        stftDuration_iterPrecent    = simCfg.stftDuration_iterPrecent;
    catch
        stftDuration_iterPrecent    = 50;
    end
    try
        snr = simCfg.snr;
    catch
        snr = inf;
    end
    try
        thetaS  = simCfg.thetaS;
    catch
        thetaS  = pi/2;
    end
    try
        targetRange_samples = simCfg.targetRange_samples;
    catch
        targetRange_samples = 32;
    end
    try
        nTheta              = simCfg.nTheta;
    catch
        nTheta              = 100;
    end
    
    propagationVelocity = 3e8;
    historyBufferSize   = 0.7;
    
    simOut.cfg.simDuration_iterations   = simDuration_iterations;
    simOut.cfg.nSensors                 = nSensors;
    simOut.cfg.r                        = r;
    simOut.cfg.compensationFreq         = compensationFreq;
    simOut.cfg.rangeErrorToLambdaRatio  = rangeErrorToLambdaRatio;
    simOut.cfg.stftDuration_iterPrecent = stftDuration_iterPrecent;
    simOut.cfg.nTheta                   = nTheta;
    simOut.cfg.targetRange_samples      = targetRange_samples;
    simOut.cfg.propagationVelocity      = propagationVelocity;
    simOut.cfg.thetaS                   = thetaS;
    simOut.cfg.historyBufferSize        = historyBufferSize;
    
    %% auxiliary
    c               = propagationVelocity;
    lambda          = c/sigFreq;
    rangeError      = lambda*rangeErrorToLambdaRatio;
    D               = lambda/2;
    f_exp           = @(f,t) exp(1i*2*pi*f*t);
    f_sig           = @(t) f_exp(sigFreq,t).*heaviside(t);
    f_noise         = @(t) rand(size(t))*10^(-snr/10);
    N               = nSensors;
    f_dTOA          = @(theta) reshape((0:(N-1))*D*cos(theta)/c,[],1);
    f_steering      = @(theta,f) fliplr(f_exp(f,f_dTOA(theta)));
    try
        simCfg.overrideFeedbackCoeffs;
        f_alphaVec  = @(thetaS,range,f,f_align) reshape(simCfg.overrideFeedbackCoeffs,[],1);
        f_betaVec   = @(thetaS,range,f,f_align) f_alphaVec(thetaS,range,f,f_align);
    catch
        f_alphaVec  = @(thetaS,range,f,f_align) ...
            reshape( ...
            (1/N) ...
            * ...
            f_steering(thetaS,f) ...
            * ...
            exp(-1i*2*pi*f_align*2*(range+rangeError)/c) ...
            ,[],1);
        f_betaVec   = @(thetaS,range,f,f_align) f_alphaVec(thetaS,range,f,f_align);
    end
    
    fSample                 = 5*sigFreq;
    tSample                 = 1/fSample;
    tPd                     = targetRange_samples*tSample;
    targetRange             = c*tPd;
    nSamplesIter            = floor(tPd/tSample);
    stftDuration_samples    = round(stftDuration_iterPrecent*nSamplesIter/100);
    firstIterStftTVec       = tSample*(0:(stftDuration_samples-1));
    stftRef                 = f_exp(-sigFreq,firstIterStftTVec);
    
    targetAngleVec          = sort(mod(thetaS + linspace(0, pi, nTheta),pi));
    simDuration_samples     = simDuration_iterations*nSamplesIter;
    hMat                    = zeros(simDuration_samples,nTheta);
    stftMat                 = zeros(simDuration_iterations,nTheta);
    
    alphaVec        = f_alphaVec(thetaS,targetRange,sigFreq,compensationFreq);
    alphaVecT       = transpose(alphaVec);
    alphaVecH       = conj(alphaVecT);
    betaVec         = f_betaVec(thetaS,targetRange,sigFreq,compensationFreq);
    betaVecT        = transpose(betaVec);
    betaVecH        = conj(betaVecT);
    sampleIdVec     = 1 : nSamplesIter;
    targetAngleId   = 0;
    for targetAngle = targetAngleVec
        targetAngleId   = targetAngleId + 1;
        dTOAVec         = f_dTOA(targetAngle);
        steerVec        = f_steering(targetAngle,sigFreq);
        feedbackSig     = zeros(simDuration_samples,1);
        arrayOutput     = zeros(simDuration_samples,1);        
        for IterId = 1 : simDuration_iterations
            iterSampleIdVec                 = (IterId-1)*nSamplesIter + sampleIdVec;
            iterTVec                        = (iterSampleIdVec-1)*tSample;
            feedbackGenerationTime          = iterTVec - 2*tPd;
            feedbackGenerationTime_samples  = round(feedbackGenerationTime/tSample);
            historySampleIdVec              = ...
                round(...
                min(feedbackGenerationTime_samples)-historyBufferSize*nSamplesIter ...
                : ...
                max(feedbackGenerationTime_samples)+historyBufferSize*nSamplesIter ...
                );
            historyTVec                     = (historySampleIdVec-1)*tSample;
            dTOAMat                         = repmat(reshape(dTOAVec,1,[]),nSamplesIter,1);
            feedbackGenerationTimeMat       = repmat(feedbackGenerationTime(:),1,N)+dTOAMat;
            curIterInput_sig                = f_sig(feedbackGenerationTimeMat);
            feedbackSig_valid               = feedbackSig(max(1,historySampleIdVec));
            curIterInput_feedback           = interp1(historyTVec(:),feedbackSig_valid(:),feedbackGenerationTimeMat,'spline');
            
            if false
                figure;
                plot(historyTVec,imag(feedbackSig_valid))
                hold on;
                for antId = 1 : size(feedbackGenerationTimeMat,2)
                    plot(feedbackGenerationTimeMat(:,antId),imag(curIterInput_feedback(:,antId)));
                end
                close all;
            end
            
            iterArrayInput                  = curIterInput_sig + r*curIterInput_feedback;
            iterArrayFeedback               = reshape(iterArrayInput*alphaVecH(:),[],1);
            feedbackSig(iterSampleIdVec)    = iterArrayFeedback(:);
            
            if false
                %% DEBUG
                figure;
                plot(abs(feedbackSig));
                close all;
            end
            
            iterArrayOutput                 = reshape(iterArrayInput*betaVecH(:),[],1) + f_noise(iterTVec(:));
            arrayOutput(iterSampleIdVec)    = iterArrayOutput(:);
            stftInput_iterH                 = iterArrayOutput(end-stftDuration_samples+1:end);
            iterHStft                       = reshape(stftRef,1,[])*stftInput_iterH(:);
            stftMat(IterId,targetAngleId)   = iterHStft;
        end
        
        hMat(:,targetAngleId)   = arrayOutput(:);
        
    end
    
    simOut.cfg.nSamplesIter = nSamplesIter;
    simOut.cfg.D            = D;
    simOut.targetAngleVec   = targetAngleVec;
    simOut.hMat             = hMat;
    simOut.stftMat          = stftMat;
end