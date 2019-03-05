function [hAbs2Rel] = H_development()
%% symbolics
syms N DU r dPhi positive;

%% generating basic terms
ad  = (r/N)*(1-exp(-1i*N*DU))/(1-exp(-1i*DU));
bd  = ad;

%% full expression transfer function
hDen    = simplify(1-ad*exp(1i*dPhi));
h       = simplify(bd/hDen);
hAbs2   = simplify(rewrite(expand(h.*conj(h)),'sincos'));

%% hOpt
hAbs2_dPhi0         = subs(hAbs2, dPhi, 0);
[lim_hAbs2,report]  = get_funcLimit(hAbs2_dPhi0,DU,0);
hAbs2Rel            = simplify(rewrite(expand(hAbs2/lim_hAbs2),'sincos'));
pretty(hAbs2Rel)
end