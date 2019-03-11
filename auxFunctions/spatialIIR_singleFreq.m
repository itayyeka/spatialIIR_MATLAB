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
nTheta              = 100;
targetRange_samples = 128;
propagationVelocity = 3e8;
thetaS              = pi/2;
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
f_sig1          = @(t) f_exp(sigFreq,t).*heaviside(t);
N               = nSensors;
f_dTOA          = @(theta) reshape(((N-1):-1:0)*D*cos(theta)/c,[],1);
try
    simCfg.overrideFeedbackCoeffs;
    f_steering      = @(theta,f) reshape(simCfg.overrideFeedbackCoeffs,[],1);
catch
    f_steering      = @(theta,f) reshape(exp(1i*2*pi*f*f_dTOA(theta)),[],1);
end
f_hCB           = @(thetaS,range,f,f_align) ...
    reshape( ...
    (1/N) ...
    * ...
    conj(f_steering(thetaS,f)) ...
    * ...
    exp(-1i*2*pi*f_align*2*(range+rangeError)/c) ...
    ,[],1);
fSample                 = 5*sigFreq;
tSample                 = 1/fSample;
tPd                     = targetRange_samples*tSample;
targetRange             = c*tPd;
nSamplesIter            = floor((2*tPd-N*D/c)/tSample);
stftDuration_samples    = round(stftDuration_iterPrecent*nSamplesIter/100);
firstIterStftTVec       = tSample*(0:(stftDuration_samples-1));
stftRef1                = f_exp(-sigFreq,firstIterStftTVec);

targetAngleVec          = linspace(0, pi, nTheta);
simDuration_samples     = simDuration_iterations*nSamplesIter;
hMat                    = zeros(simDuration_samples,nTheta);
stftMat                 = zeros(simDuration_iterations,nTheta);

hCB             = f_hCB(thetaS,targetRange,sigFreq,compensationFreq);
hCBT            = transpose(hCB);
sampleIdVec     = 1 : nSamplesIter;
targetAngleId   = 0;
for targetAngle = targetAngleVec
    targetAngleId   = targetAngleId + 1;
    dTOAVec         = f_dTOA(targetAngle);
    h               = [];
    
    for IterId = 1 : simDuration_iterations
        iterSampleIdVec             = (IterId-1)*nSamplesIter + sampleIdVec;
        iterTVec                    = (iterSampleIdVec-1)*tSample;
        historyStartIterId          = IterId-(1+historyIterNum);
        historyEndIterId            = IterId-1;
        historySampleIdVec          = ((1+historyStartIterId*nSamplesIter) : historyEndIterId*nSamplesIter);
        historyTVec                 = (historySampleIdVec-1)*tSample;
        feedbackGenerationTime      = iterTVec - 2*tPd;
        feedbackGenerationTimeMat   = repmat(feedbackGenerationTime(:),1,N)-repmat(reshape(dTOAVec,1,[]),nSamplesIter,1);
        
        curIterInput_sig                = f_sig1(feedbackGenerationTimeMat);
        curIterInput_feedback           = f_resample(historyTVec(:),h(:),feedbackGenerationTimeMat);
        
        iterH = reshape( ...
            ( ...
            curIterInput_sig ...
            + ...
            r*curIterInput_feedback ...
            ) ...
            * ...
            hCBT(:) ...
            ,[],1);
        
        h = ...
            [...
            h(:) ...
            ; ...
            iterH(:) ...
            ];
        
        if false
            %% DEBUG
            figure; plot(real(h1));
            close all;
        end
        
        stftInput_iterH                 = iterH(end-stftDuration_samples+1:end);
        iterHStft                       = reshape(stftRef1,1,[])*stftInput_iterH(:);
        stftMat(IterId,targetAngleId)   = iterHStft;
    end
    
    hMat(:,targetAngleId)   = h(:);
    
end

simOut.targetAngleVec   = targetAngleVec;
simOut.hMat             = hMat;
simOut.stftMat          = stftMat;
end