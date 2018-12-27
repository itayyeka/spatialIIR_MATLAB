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
defautCfg.backoffFactor_min             = 0;
defautCfg.backoffFactor_max             = 1;
defautCfg.backoffFactor_nValues         = 10;
defautCfg.rErr_min                      = -20;
defautCfg.rErr_max                      = 20;
defautCfg.rErr_nValues                  = 30;
defautCfg.thetaSim_azimuthalWidth       = pi/4;
defautCfg.thetaSim_nPoints              = 5;
defautCfg.polesThetaVec                 = [pi/3 -pi/3];
defautCfg.RangeVal                      = 1000;
defautCfg.DVal                          = 0.01;
defautCfg.cVal                          = 3e8;
defautCfg.fSample                       = 100e6;
defautCfg.syncSigCfg.duration_SAMPLES   = 1024;
defautCfg.syncSigCfg.baseFreq_relative  = 0.1;
defautCfg.syncSigCfg.bandwidth_relative = 0.1;
defautCfg.integralNPoints               = 2048;

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
end

simOutput = [];
structToVariables(crlbSymbolics);

%% clculate the spectral content of s(t)
sData           = generateSignalData(syncSigCfg);
sSpectrum       = fftshift(fft(sData));
sSpectrumFreq   = linspace(-0.5,0.5,syncSigCfg.duration_SAMPLES);
%% prepare for numerical integration
integralFreqVec             = linspace(-0.5,0.5,integralNPoints);
sSpectrumInterpulated       = interp1(sSpectrumFreq, sSpectrum, integralFreqVec, 'spline');
sSpectrumInterpulated_abs2  = sSpectrumInterpulated.*conj(sSpectrumInterpulated);
dOmega                      = 2*pi/integralNPoints;

%% simulating
if true
    %% creating the ideal beta values
    syms                x;
    poles               = exp(1i*omega*D*cos(polesThetaVec)/c);
    polynomBeta         = expand(prod(x-poles));
    betaV_ideal         = coeffs(polynomBeta,x);
    phaseAlign          = exp(1i*omega*tau);
    betaOpt_sym         = sym(zeros(size(betaV_ideal)));
    betaOpt_sym(1)      = (1-betaV_ideal(1))*phaseAlign;
    betaOpt_sym(2:end)  = -betaV_ideal(2:end)*phaseAlign;
    
    %% setting constant parameter values
    alphaValues     = zeros(size(alphaV));
    alphaValues(1)  = 1;
    thetaVal        = polesThetaVec(1);
    
    %% prepare parameters grid
    backoffFactorVec                        = linspace(backoffFactor_min    ,backoffFactor_max  ,backoffFactor_nValues);
    rErrVec                                 = linspace(rErr_min             ,rErr_max           ,rErr_nValues);
    thetaVec                                = thetaVal+linspace(-0.5,0.5,thetaSim_nPoints)*thetaSim_azimuthalWidth;
    [backoffFactorMat, rErrMat, thetaMat]   = meshgrid(backoffFactorVec, rErrVec, thetaVec);
    crlb_theta_resultMat                           = zeros(size(backoffFactorMat));
    dynamicPrmSet                           = [];
    
    %% generating symbolic & values vec for easy substitution
    symVarVec_noBeta_noOmega    = [D    ;c      ;targetRange    ;theta      ;alphaV(:)  ];
    symVarVec_noBeta            = [D    ;c      ;targetRange    ;theta      ;alphaV(:)      ;omega              ];
    symVarVec_noOmega           = [symVarVec_noBeta_noOmega(:) ; betaV(:)];
    omegaVal_betaCalc           = 2*pi*fSample*syncSigCfg.baseFreq_relative;
    symVarValues_forBetaCalc    = [DVal ;cVal   ;RangeVal       ;thetaVal   ;alphaValues(:) ;omegaVal_betaCalc  ];
    betaValues                  = subs(betaOpt_sym, symVarVec_noBeta, symVarValues_forBetaCalc);
    
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
    
    %% calculating CRLB values
    paramSetSize    = size(crlb_theta_resultMat);
    omegaSymMat     = repmat(omega, size(crlb_theta_resultMat));
    paramSet_CELL   = cell(paramSetSize);
    for paramSetID  = 1:numel(crlb_theta_resultMat)
        paramSet_CELL{paramSetID}   = [DVal     ;cVal    ;RangeVal + rErrMat(paramSetID)    ;thetaMat(paramSetID)   ;alphaValues(:) ;betaValues(:)];
    end
    symInptPartSpectum      = cell(paramSetSize);
    integralResultVec       = cell(paramSetSize);
    integralResult          = cell(paramSetSize);
    integralFreqVec_CELL    = repmat({integralFreqVec}, paramSetSize);
    
    parfor paramSetID = 1:numel(paramSet_CELL)        
        %% symbolic parts of integrands
        if true
            %% J_theta_theta
            curSymIntPart   = subs(J_theta_theta_symIntPart,symVarVec_noOmega,paramSet_CELL{paramSetID});
            
            symInptPartSpectum{paramSetID}  = subs(curSymIntPart,omegaSymMat(paramSetID),integralFreqVec_CELL{paramSetID});
            integralResultVec{paramSetID}   = symInptPartSpectum{paramSetID}.*sSpectrumInterpulated_abs2;
            integralResult{paramSetID}      = sum(integralResultVec{paramSetID});
            
            J_theta_theta   = integralResult{paramSetID};
            %% J_tau_tau
            curSymIntPart   = subs(J_tau_tau_symIntPart,symVarVec_noOmega,paramSet_CELL{paramSetID});
            
            symInptPartSpectum{paramSetID}  = subs(curSymIntPart,omegaSymMat(paramSetID),integralFreqVec_CELL{paramSetID});
            integralResultVec{paramSetID}   = symInptPartSpectum{paramSetID}.*sSpectrumInterpulated_abs2;
            integralResult{paramSetID}      = sum(integralResultVec{paramSetID});
            
            J_tau_tau       = integralResult{paramSetID};
            %% J_theta_tau
            curSymIntPart       = subs(J_theta_tau_symIntPart,symVarVec_noOmega,paramSet_CELL{paramSetID});
            
            symInptPartSpectum{paramSetID}  = subs(curSymIntPart,omegaSymMat(paramSetID),integralFreqVec_CELL{paramSetID});
            integralResultVec{paramSetID}   = symInptPartSpectum{paramSetID}.*sSpectrumInterpulated_abs2;
            integralResult{paramSetID}      = sum(integralResultVec{paramSetID});
            
            J_theta_tau         = integralResult{paramSetID};
            J_tau_theta         = J_theta_tau;
        end
        
        J = [J_theta_theta J_theta_tau ; J_tau_tau J_tau_theta];
        
        Jinv = inv(J);
        
        %% collect result
        crlb_theta_resultMat(paramSetID)   = Jinv(1,1);
    end
    
    crlb_theta_resultMat   = eval('crlbResultMat;');
    CRLB_theta_norm = crlb_theta_resultMat / max(abs(crlb_theta_resultMat(:)));
end
figure;surf(backoffFactorMat, rErrMat, db(abs(CRLB_theta_norm)));
end