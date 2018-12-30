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
integralFreqVec             = linspace(-0.5,0.5,integralNPoints);
dOmega                      = 2*pi/integralNPoints;

%% auxiliary symbolics
betaValues                  = subs(betaOpt_sym, ...
    [ D,    c,      omega,                                  targetRange ], ...
    [ DVal, cVal,   fSample*syncSigCfg.baseFreq_relative,   RangeVal    ] ...
    );
alphaValues                 = zeros(size(alphaV));
alphaValues(1)              = 1;

%% calculate FIM related symbolic expressions
if true
    %% J_theta_theta_symIntPart
    term1           = (1-g)^2;
    term2           = ...
        (alphaT*(A*d+B*betaV*exp(-1i*omega*tau))) ...
        / ...
        term1;
    
    J_theta_theta_symIntPart    = term2*conj(term2);
    %% J_tau_tau_symIntPart
    term3           = ...
        (alphaT*d) ...
        / ...
        term1;
    
    J_tau_tau_symIntPart    = (omega^2)*term3*conj(term3);
    %% J_theta_tau_symIntPart
    term4           = ...
        (1i*omega*alphaT*(A*d+B*betaV*exp(-1i*omega*tau))*alphaH*conj(d));
    
    J_theta_tau_symIntPart  = real(term4/term1);
end

%% prepare parameters grid
backoffFactorVec    = linspace(backoffFactor_min    ,backoffFactor_max  ,backoffFactor_nValues);
rErr_min            = 0;
baseferqLambda      = cVal/(syncSigCfg.baseFreq_relative*fSample);
rErr_max            = baseferqLambda;
rErrVec             = linspace(rErr_min             ,rErr_max           ,rErr_nValues);
thetaVec            = polesThetaVec(1)+linspace(-0.5,0.5,thetaSim_nPoints)*thetaSim_azimuthalWidth;
bandwidthVec        = fliplr(linspace(syncSigCfg.bandwidth_relative_max,0,syncSigCfg.bandwidth_relative_nValues));
simOutData_CELL     = cell(backoffFactor_nValues * rErr_nValues * thetaSim_nPoints * syncSigCfg.bandwidth_relative_nValues, 1);

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
if true
    %% stage 1
    parfor outputID = 1 : numel(simOutData_CELL)
        disp(['tescase #' num2str(outputID)]);
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
        sSpectrumFreq                   = linspace(-0.5,0.5,generateSignalData_cfgIn.duration_SAMPLES);
        
        %% prepare for numerical integration
        sSpectrumInterpulated       = interp1(sSpectrumFreq, sSpectrum, integralFreqVec, 'spline');
        sSpectrumInterpulated_abs2  = sSpectrumInterpulated.*conj(sSpectrumInterpulated);
        
        if useSimplifiedSignal
            sigFreqSupport_relative     = ...
                syncSigCfg.baseFreq_relative ...
                +...
                [-0.5 0.5] * syncSigCfg.bandwidth_relative_max;
            sSpectrumInterpulated_abs2  = ...
                double(abs(integralFreqVec)>=sigFreqSupport_relative(1)) ...
                .* ...
                double(abs(integralFreqVec)<=sigFreqSupport_relative(2)) ...
                ;
            
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
                curSymIntPart   = subs(J_theta_theta_symIntPart,J_symVarVec,J_symVarValues);
                
                symInptPartSpectum  = simplify(subs(curSymIntPart,omega,integralFreqVec));
                integralResultVec   = symInptPartSpectum.*sSpectrumInterpulated_abs2;
                integralResult      = sum(integralResultVec);
                
                J_theta_theta   = integralResult;
                %% J_tau_tau
                curSymIntPart   = subs(J_tau_tau_symIntPart,J_symVarVec,J_symVarValues);
                
                symInptPartSpectum  = subs(curSymIntPart,omega,integralFreqVec);
                integralResultVec   = symInptPartSpectum.*sSpectrumInterpulated_abs2;
                integralResult      = sum(integralResultVec);
                
                J_tau_tau       = integralResult;
                %% J_theta_tau
                curSymIntPart       = subs(J_theta_tau_symIntPart,J_symVarVec,J_symVarValues);
                
                symInptPartSpectum  = subs(curSymIntPart,omega,integralFreqVec);
                integralResultVec   = symInptPartSpectum.*sSpectrumInterpulated_abs2;
                integralResult      = sum(integralResultVec);
                
                J_theta_tau         = integralResult;
                J_tau_theta         = J_theta_tau;
            end
            %% form the information matrix
            J       = dOmega*[J_theta_theta J_theta_tau ; J_tau_tau J_tau_theta];
        end
        
        simOutData_CELL{outputID}.Jsym  = J;
    end
end

warning('on');

extractedSimOutData_CELL    = cell(size(simOutData_CELL));
parfor outputID = 1 : numel(simOutData_CELL)
    disp(['extracting results from tescase #' num2str(outputID)]);
    extractedSimOutData_CELL{outputID}  = f_extractDataFromSimOut(simOutData_CELL{outputID});
end
simOutput.simOutData_CELL = extractedSimOutData_CELL;

end

function [CELL_out] = f_extractDataFromSimOut(CELL_in)
CELL_out.paramSet   = CELL_in.paramSet;
CELL_out.J          = eval(CELL_in.Jsym);
CELL_out.Jinv       = inv(CELL_out.J);
end