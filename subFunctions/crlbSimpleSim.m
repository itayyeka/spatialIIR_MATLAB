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

%% creating the beta values
syms                        x;
polesThetaVec               = [pi/3 -pi/3]; 
poles                       = exp(1i*omega*D*cos(polesThetaVec)/c); 
polynomBeta                 = expand(prod(x-poles));
betaV_ideal                 = coeffs(polynomBeta,x);
phaseAlign                  = exp(1i*omega*tau);
betaV_values_sym            = sym(zeros(size(betaV_ideal)));
betaV_values_sym(1)         = (1-betaV_ideal(1))*phaseAlign;
betaV_values_sym(2:end)     = -betaV_ideal(2:end)*phaseAlign;

%% simulating
alphaValues     = zeros(size(alphaV));
alphaValues(1)  = 1;
RangeVal        = 1000;
omegaVal        = 2*pi*1e6;
DVal            = 0.01;
cVal            = 3e8;
thetaVal        = polesThetaVec(1);

symVarVec       = [D    ; omega     ; c     ; alphaV(:)      ; Range     ; theta      ];
symVarValuesVec = [DVal ; omegaVal  ; cVal  ; alphaValues(:) ; RangeVal  ; thetaVal   ];

CRLB_theta_kappa_sym    = subs(CRLB_theta,              symVarVec, symVarValuesVec);
betaV_values            = reshape(eval(subs(betaV_values_sym,   symVarVec, symVarValuesVec)),size(betaV));

pretty(CRLB_theta_kappa_sym)

CRLB_theta_kappa    = [];
for kappa = 0.1 : 0.01 : 1
    CRLB_theta_kappa(end+1)    = eval(subs(CRLB_theta_kappa_sym(end),  betaV, kappa*betaV_values));
end
CRLB_theta_norm = CRLB_theta_kappa / max(abs(CRLB_theta_kappa));
figure;plot(abs(CRLB_theta_norm(:)));
end