function [simOutput] = crlbSimpleSim(crlbSymbolics,cfgIn)
%{
Assuming:
s(w)    = delta(w-w0)
alphaV  = [1 0 ... 0]^T

and setting:
c       = 3e8;
Range   = 1000;
D       = 0.01;
theta   = pi/4; (this is where the poles are set)
omega   = 2*pi*1e6;
%}
if ~isstruct(cfgIn)
    cfgIn = struct();
end
%% configure
defautCfg.backoffFactor_min                         = 0.7;
defautCfg.backoffFactor_max                         = 0.8;
defautCfg.backoffFactor_nValues                     = 10;
defautCfg.rErr_nValues                              = 10;
defautCfg.thetaSim_azimuthalWidth                   = 0;
defautCfg.thetaSim_nPoints                          = 2;
defautCfg.polesThetaVec                             = [pi/3 -pi/3];
defautCfg.RangeVal                                  = 1000;
defautCfg.DVal                                      = 0.01;
defautCfg.cVal                                      = 3e8;
defautCfg.fSample                                   = 100e6;
defautCfg.syncSigCfg.duration_SAMPLES               = 1024;
defautCfg.syncSigCfg.baseFreq_relative              = 0.1;
defautCfg.syncSigCfg.bandwidth_relative_max         = 0.1;
defautCfg.syncSigCfg.bandwidth_relative_nValues     = 10;
defautCfg.integralNPoints                           = 128;
defautCfg.useSimplifiedSignal                       = 1;

cfgFields = fieldnames(defautCfg);

for cfgFieldID = 1 : numel(cfgFields)
    curCfgField     = cfgFields{cfgFieldID};
    try
        cmdString   = [curCfgField '=cfgIn.(''' curCfgField ''');'];
        eval(cmdString);
    catch
        cmdString   = [curCfgField '=defautCfg.(''' curCfgField ''');'];
        eval(cmdString);
    end
    finalCfg.(curCfgField) = eval([curCfgField ';']);
end

simOutput       = [];
simOutput.cfg   = finalCfg;
structToVariables(crlbSymbolics);

%% creating the ideal beta values
syms                x;
poles               = exp(1i*omega*D*cos(polesThetaVec)/c);
polynomBeta         = expand(prod(x-poles));
betaV_ideal         = coeffs(polynomBeta,x);
phaseAlign          = exp(1i*omega*tau);
betaOpt_sym         = sym(zeros(size(betaV_ideal)));
betaOpt_sym(1)      = (1-betaV_ideal(1))*phaseAlign;
betaOpt_sym(2:end)  = -betaV_ideal(2:end)*phaseAlign;

%% integral auxiliary veriables
integralFreqVec_relative    = linspace(-0.5,0.5,integralNPoints);
integralOmegaVec            = 2*pi*fSample*integralFreqVec_relative;
dOmega                      = 2*pi/integralNPoints;

%% auxiliary symbolics
betaValues                  = subs(betaOpt_sym, ...
    [ D,    c,      omega,                                      targetRange ], ...
    [ DVal, cVal,   2*pi*fSample*syncSigCfg.baseFreq_relative,  RangeVal    ] ...
    );
alphaValues                 = zeros(size(alphaV));
alphaValues(1)              = 1;

%% calculate FIM related symbolic expressions
if true
    J_symVarVec     =   [ D;        alphaV(:);      betaV(:);                           c;      targetRange;        theta;  omega];
    %% J_theta_theta_symIntPart
    term1           = (1-g)^2;
    term2           = ...
        (alphaT*(A*d+B*betaV*exp(-1i*omega*tau))) ...
        / ...
        term1;
    
    J_theta_theta_symIntPart    = term2*conj(term2);
    f_J_theta_theta             = matlabFunction(simplify(J_theta_theta_symIntPart), 'Vars', {[J_symVarVec]});
    %% J_tau_tau_symIntPart
    term3           = ...
        (alphaT*d) ...
        / ...
        term1;
    
    J_tau_tau_symIntPart    = (omega^2)*term3*conj(term3);
    f_J_tau_tau             = matlabFunction(simplify(J_tau_tau_symIntPart), 'Vars', {[J_symVarVec]});
    %% J_theta_tau_symIntPart
    term4           = ...
        (1i*omega*alphaT*(A*d+B*betaV*exp(-1i*omega*tau))*alphaH*conj(d));
    
    J_theta_tau_symIntPart  = real(term4/term1);
    f_J_theta_tau           = matlabFunction(simplify(J_theta_tau_symIntPart), 'Vars', {[J_symVarVec]});
end

%% prepare parameters grid
backoffFactorVec    = linspace(backoffFactor_min    ,backoffFactor_max  ,backoffFactor_nValues);
rErr_min            = 0;
baseferqLambda      = cVal/(syncSigCfg.baseFreq_relative*fSample);
rErr_max            = 2*baseferqLambda;
rErrVec             = linspace(rErr_min             ,rErr_max           ,rErr_nValues);
thetaVec            = polesThetaVec(1)+linspace(-0.5,0.5,thetaSim_nPoints)*thetaSim_azimuthalWidth;
bandwidthVec        = fliplr(linspace(syncSigCfg.bandwidth_relative_max,0,syncSigCfg.bandwidth_relative_nValues));
simOutData_CELL     = cell(backoffFactor_nValues * rErr_nValues * thetaSim_nPoints * syncSigCfg.bandwidth_relative_nValues, 1);

simOutput.cfg.prmVec.bandwidthVec       = bandwidthVec;
simOutput.cfg.prmVec.backoffFactorVec   = backoffFactorVec;
simOutput.cfg.prmVec.rErrVec            = rErrVec;
simOutput.cfg.prmVec.thetaVec           = thetaVec;

noFeedbackSimulationFlag        = 0;
if isequal(bandwidthVec,0)
    noFeedbackSimulationFlag        = 1;
    polesThetaVec
    noFeedbackSimAlpha_symVarVec    = [D;       c;      omega;                                      theta];
    noFeedbackSimAlpha_symVarValues = [DVal;    cVal;   2*pi*fSample*syncSigCfg.baseFreq_relative;  polesThetaVec(1)];
    alphaValues_preNorm             = conj(eval(subs(d,noFeedbackSimAlpha_symVarVec,noFeedbackSimAlpha_symVarValues)));    
    alphaValues                     = alphaValues_preNorm/norm(alphaValues_preNorm);
end

outputID = 1;
for bandwidthID = 1 : syncSigCfg.bandwidth_relative_nValues
    for backoffFactorID = 1 : backoffFactor_nValues
        for rErrID = 1 : rErr_nValues
            for thetaID = 1 : thetaSim_nPoints
                paramSet.bandwidth                  = bandwidthVec(bandwidthID);
                paramSet.backoffFactor              = backoffFactorVec(backoffFactorID);
                paramSet.rErr                       = rErrVec(rErrID);
                paramSet.theta                      = thetaVec(thetaID);
                simOutData_CELL{outputID}.paramSet  = paramSet;
                outputID                            = outputID + 1;
            end
        end
    end
end

%% simulating
if true
    %% workspace handling related assigning
    syncSigCfg              = syncSigCfg;
    D                       = D;
    alphaV                  = alphaV;
    betaV                   = betaV;
    c                       = c;
    targetRange             = targetRange;
    theta                   = theta;
    DVal                    = DVal;
    alphaValues             = alphaValues;
    betaValues              = betaValues;
    cVal                    = cVal;
    RangeVal                = RangeVal;
    omega                   = omega;
    useSimplifiedSignal     = useSimplifiedSignal;
end

warning('off');

bandiwdthEffectDEBUG    = 0;
if bandiwdthEffectDEBUG %% DEBUG
    %% DEBUG
    figure;
    hold on;
end

parfor outputID = 1 : numel(simOutData_CELL)
    curParamSet         = simOutData_CELL{outputID}.paramSet;
    curBandwidth        = curParamSet.bandwidth;
    curBackoffFactor    = curParamSet.backoffFactor;
    curRErr             = curParamSet.rErr;
    curTheta            = curParamSet.theta;
    
    %% clculate the spectral content of s(t)
    generateSignalData_cfgIn                        = syncSigCfg;
    generateSignalData_cfgIn.bandwidth_relative     = curBandwidth;
    
    sData                           = generateSignalData(generateSignalData_cfgIn);
    sSpectrum                       = fftshift(fft(sData));
    sSpectrumFreq_relative          = linspace(-0.5,0.5,generateSignalData_cfgIn.duration_SAMPLES);
    
    %% prepare for numerical integration
    sSpectrumInterpulated       = interp1(sSpectrumFreq_relative,   sSpectrum, integralFreqVec_relative, 'spline');
    sSpectrumInterpulated_abs2  = sSpectrumInterpulated.*conj(sSpectrumInterpulated);
    
    if useSimplifiedSignal
        sigFreqSupport_relative     = ...
            syncSigCfg.baseFreq_relative ...
            +...
            [-0.5 0.5] * curBandwidth;
        sSpectrumInterpulated_abs2  = ...
            double(abs(integralFreqVec_relative)>=sigFreqSupport_relative(1)) ...
            .* ...
            double(abs(integralFreqVec_relative)<=sigFreqSupport_relative(2)) ...
            ;
        
        if ~any(sSpectrumInterpulated_abs2)
            %% for zero bnadwidth scenarios
            [freqVal,   freqID      ] = min(abs(syncSigCfg.baseFreq_relative-integralFreqVec_relative));
            [~,         negFreqID   ] = min(abs(-freqVal                    -integralFreqVec_relative));
            sSpectrumInterpulated_abs2([negFreqID, freqID])     = 1;
        end
        
        if false
            %% DEBUG
            figure; plot(sSpectrumInterpulated_abs2);
        end
    end
    
    %% calculating CRLB values
    if true
        %% symbolic parts of integrands
        if true
            J_symVarVec     =   [ D;        alphaV(:);      betaV(:);                           c;      targetRange;        theta   ];
            J_symVarValues  =   [ DVal;     alphaValues(:); curBackoffFactor*betaValues(:);     cVal;   RangeVal + curRErr; curTheta];
            %% J_theta_theta
            curSymIntPart       = subs(J_theta_theta_symIntPart,J_symVarVec,J_symVarValues);
            
            symInptPartSpectum  = subs(curSymIntPart,omega,integralOmegaVec);           
            %symInptPartSpectum  = cellfun(@(curOmega) f_J_theta_theta([J_symVarValues(:) ; curOmega]), num2cell(integralOmegaVec));
            integralResultVec   = symInptPartSpectum.*sSpectrumInterpulated_abs2;
            integralResult      = sum(integralResultVec);
            
            if bandiwdthEffectDEBUG %% DEBUG
                %% DEBUG
                plot([real(eval(integralResultVec))]);
            end
            
            J_theta_theta   = integralResult;
            %% J_tau_tau
            curSymIntPart   = subs(J_tau_tau_symIntPart,J_symVarVec,J_symVarValues);
            
            symInptPartSpectum  = subs(curSymIntPart,omega,integralOmegaVec);
            integralResultVec   = symInptPartSpectum.*sSpectrumInterpulated_abs2;
            integralResult      = sum(integralResultVec);
            
            J_tau_tau       = integralResult;
            
            %% J_theta_tau
            curSymIntPart       = subs(J_theta_tau_symIntPart,J_symVarVec,J_symVarValues);
            
            symInptPartSpectum  = subs(curSymIntPart,omega,integralOmegaVec);
            integralResultVec   = symInptPartSpectum.*sSpectrumInterpulated_abs2;
            integralResult      = sum(integralResultVec);
            
            J_theta_tau         = integralResult;
            J_tau_theta         = J_theta_tau;
        end
        %% form the information matrix
        J       = dOmega*[J_theta_theta J_theta_tau ; J_tau_tau J_tau_theta];
    end    
    simOutData_CELL{outputID}.Jsym  = J;
    disp(['finished calculation of tescase #' num2str(outputID)]);
end
disp(['finished calculation of all test cases']);
warning('on');
disp(['Initiating a container for the sim results.']);
extractedSimOutData_CELL    = cell(size(simOutData_CELL));
parfor outputID = 1 : numel(simOutData_CELL)
    extractedSimOutData_CELL{outputID}  = f_extractDataFromSimOut(simOutData_CELL{outputID});
    disp(['finished data extraction from tescase #' num2str(outputID)]);
end
disp(['extracted data from all test cases']);
simOutput.simOutData_CELL = extractedSimOutData_CELL;

end

function [CELL_out] = f_extractDataFromSimOut(CELL_in)
CELL_out.paramSet   = CELL_in.paramSet;
CELL_out.J          = eval(CELL_in.Jsym);
CELL_out.Jinv       = inv(CELL_out.J);
end