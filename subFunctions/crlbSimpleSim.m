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
%% configure
backoffFactor_min       = 0;
backoffFactor_max       = 1;
backoffFactor_nValues   = 10;
rErr_min                = -20;
rErr_max                = 20;
rErr_nValues            = 30;

%%

simOutput = [];
structToVariables(crlbSymbolics);

term1           = (1-g)^2;
term2           = ...
    (alphaT*(A*d+B*betaV*exp(-1i*omega*tau))) ...
    / ...
    term1;

J_theta_theta   = term2*conj(term2);

term3           = ...
    (alphaT*d) ...
    / ...
    term1;

J_tau_tau       = (omega^2)*term3*conj(term3);

term4           = ...
    (1i*omega*alphaT*(A*d+B*betaV*exp(-1i*omega*tau))*alphaH*conj(d));

J_theta_tau     = real(term4/term1);
J_tau_theta     = J_theta_tau;

J = [J_theta_theta J_theta_tau ; J_tau_tau J_tau_theta];

Jinv = inv(J);

CRLB_theta = Jinv(1,1);

%% simulating
if true
    %% creating the ideal beta values
    syms                x;
    polesThetaVec       = [pi/3 -pi/3];
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
    RangeVal        = 1000;
    omegaVal        = 2*pi*1e6;
    DVal            = 0.01;
    cVal            = 3e8;
    thetaVal        = polesThetaVec(1);
    
    %% generating symbolic & values vec for easy substitution
    symVarVec               = [D    ; omega     ; c     ; alphaV(:)      ; targetRange      ; theta      ];
    symVarValuesVecIdeal    = [DVal ; omegaVal  ; cVal  ; alphaValues(:) ; RangeVal         ; thetaVal   ];
    syms                    rErr;
    symVarValuesVec         = [DVal ; omegaVal  ; cVal  ; alphaValues(:) ; RangeVal + rErr  ; thetaVal   ];
        
    %% calculating CRLB values
    betaOptIdeal        = eval(subs(betaOpt_sym,  symVarVec,  symVarValuesVecIdeal));
    CRLB_theta_sym      = subs(CRLB_theta,  symVarVec,  symVarValuesVec);    
    backoffFactorVec    = linspace(backoffFactor_min    ,backoffFactor_max  ,backoffFactor_nValues);
    rErrVec             = linspace(rErr_min             ,rErr_max           ,rErr_nValues);
    
    [backoffFactorMat, rErrMat] = meshgrid(backoffFactorVec, rErrVec);
    crlbResultMat               = zeros(size(backoffFactorMat));
    dynamicPrmSet               = [betaV(:) ; rErr];
    
    for paramSetID = 1:numel(crlbResultMat)
        crlbResultMat(paramSetID)   = ...
            subs(...
            CRLB_theta_sym, ...
            dynamicPrmSet, ...
            [backoffFactorMat(paramSetID)*betaOptIdeal(:) ; rErrMat(paramSetID)] ...
            ) ...
            ;
    end
    
    crlbResultMat   = eval('crlbResultMat;');
    CRLB_theta_norm = crlbResultMat / max(abs(crlbResultMat(:)));
end
figure;surf(backoffFactorMat, rErrMat, db(abs(CRLB_theta_norm)));
end