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
historyIterNum      = 3;

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
simOut.cfg.historyIterNum           = historyIterNum;

%% auxiliary
c               = propagationVelocity;
lambda          = c/sigFreq;
rangeError      = lambda*rangeErrorToLambdaRatio;
D               = lambda/2;
f_exp           = @(f,t) exp(-1i*2*pi*f*t);
f_sig           = @(t) f_exp(sigFreq,t).*heaviside(t);
f_noise         = @(t) rand(size(t))*10^(-snr/10);
N               = nSensors;
f_dTOA          = @(theta) reshape(((N-1):-1:0)*D*cos(theta)/c,[],1);
f_steering      = @(theta,f) reshape(exp(1i*2*pi*f*f_dTOA(theta)),[],1);
try
    simCfg.overrideFeedbackCoeffs;
    f_alphaVec  = @(thetaS,range,f,f_align) reshape(simCfg.overrideFeedbackCoeffs,[],1);
    f_betaVec   = @(thetaS,range,f,f_align) f_alphaVec(thetaS,range,f,f_align);
catch
    f_alphaVec  = @(thetaS,range,f,f_align) ...
        reshape( ...
        (1/N) ...
        * ...
        conj(f_steering(thetaS,f)) ...
        * ...
        exp(-1i*2*pi*f_align*2*(range+rangeError)/c) ...
        ,[],1);
    f_betaVec   = @(thetaS,range,f,f_align) f_alphaVec(thetaS,range,f,f_align);
end

fSample                 = 5*sigFreq;
tSample                 = 1/fSample;
tPd                     = targetRange_samples*tSample;
targetRange             = c*tPd;
nSamplesIter            = floor((2*tPd-N*D/c)/tSample);
stftDuration_samples    = round(stftDuration_iterPrecent*nSamplesIter/100);
firstIterStftTVec       = tSample*(0:(stftDuration_samples-1));
stftRef                 = f_exp(-sigFreq,firstIterStftTVec);

targetAngleVec          = linspace(0, pi, nTheta);
simDuration_samples     = simDuration_iterations*nSamplesIter;
hMat                    = zeros(simDuration_samples,nTheta);
stftMat                 = zeros(simDuration_iterations,nTheta);

alphaVec        = f_alphaVec(thetaS,targetRange,sigFreq,compensationFreq);
alphaVecT       = transpose(alphaVec);
betaVec         = f_betaVec(thetaS,targetRange,sigFreq,compensationFreq);
betaVecT        = transpose(betaVec);
sampleIdVec     = 1 : nSamplesIter;
targetAngleId   = 0;
for targetAngle = targetAngleVec
    targetAngleId   = targetAngleId + 1;
    dTOAVec         = f_dTOA(targetAngle);
    feedbackSig     = zeros(simDuration_samples,1);
    arrayOutput     = zeros(simDuration_samples,1);
    
    for IterId = 1 : simDuration_iterations
        iterSampleIdVec             = (IterId-1)*nSamplesIter + sampleIdVec;
        iterTVec                    = (iterSampleIdVec-1)*tSample;
        historyStartIterId          = IterId-(1+historyIterNum);
        historyEndIterId            = IterId-1;
        historySampleIdVec          = ((1+historyStartIterId*nSamplesIter) : historyEndIterId*nSamplesIter);
        historyTVec                 = (historySampleIdVec-1)*tSample;
        historySampleIdVec_valid    = historySampleIdVec(historySampleIdVec>=1);
        historyTVec_valid           = historyTVec(historySampleIdVec>=0);
        feedbackGenerationTime      = iterTVec - 2*tPd;
        feedbackGenerationTimeMat   = repmat(feedbackGenerationTime(:),1,N)-repmat(reshape(dTOAVec,1,[]),nSamplesIter,1);
        
        if isinf(snr)
            curIterInput_sig    = f_sig(feedbackGenerationTimeMat);
        else
            curIterInput_sig    = ...
                f_sig(feedbackGenerationTimeMat) ...
                ...+ ...
                ...f_noise(feedbackGenerationTimeMat) ...
                ;
        end
        
        curIterInput_feedback           = ...
            f_resample(historyTVec_valid,feedbackSig(historySampleIdVec_valid),feedbackGenerationTimeMat) ...
            ...+ ...
            ...f_noise(feedbackGenerationTimeMat)
            ;
        iterArrayInput                  = curIterInput_sig + r*curIterInput_feedback;
        iterArrayFeedback               = reshape(iterArrayInput*alphaVecT(:),[],1);
        feedbackSig(iterSampleIdVec)    = iterArrayFeedback(:);
        
        if false
            %% DEBUG
            figure;plot(real(feedbackSig));
            close all;            
        end
        
        iterArrayOutput                 = reshape(iterArrayInput*betaVecT(:),[],1) + f_noise(iterTVec(:));
        arrayOutput(iterSampleIdVec)    = iterArrayOutput(:);
        stftInput_iterH                 = iterArrayOutput(end-stftDuration_samples+1:end);
        iterHStft                       = reshape(stftRef,1,[])*stftInput_iterH(:);
        stftMat(IterId,targetAngleId)   = iterHStft;
    end
    
    hMat(:,targetAngleId)   = arrayOutput(:);
    
end

simOut.targetAngleVec   = targetAngleVec;
simOut.hMat             = hMat;
simOut.stftMat          = stftMat;
end