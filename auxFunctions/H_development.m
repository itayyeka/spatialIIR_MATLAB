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

%% analytical
DUVec       = linspace(-pi+1e-6,pi-1e-6,1000);
hSteerErr   = subs(h,dPhi,0);
if true
    %% HPBW
    curExpr     = diff(diff(subs(hAbs2Rel,dPhi,0),DU),DU);
    nTaylor     = 4;    
    taylorPoly  = 1;
    nVal        = 5;
    rVal        = 0.7;
    syms X;
    diffExpr    = cell(nTaylor,1);
    if true
        for diffID = 2 : 2 : 2*nTaylor
            diffExpr{diffID/2}  = get_funcLimit(curExpr,DU,0);
            curExpr             = diff(diff(curExpr,DU),DU);
            taylorPoly          = taylorPoly + diffExpr{diffID/2}*(X^diffID)/factorial(diffID);
        end
    end    
    taylorPoly_sim  = abs(eval(subs(taylorPoly,{N,r,X}, {nVal,rVal,DUVec})));
    hSteerErr_sim   = abs(eval(subs(hAbs2Rel,{N,r,DU,dPhi}, {nVal,rVal,DUVec,0})));
    figure;
    hold on;
    plot(DUVec,taylorPoly_sim);
    plot(DUVec,hSteerErr_sim,'--'); 
    ylim([0 1]);
    %% sidelobes
    if false        
        [~,hSteerErr_den]   = numden(hSteerErr);
        curExpr             = (1/hSteerErr_den)*conj(1/hSteerErr_den);
        curExpr_sim         = abs(eval(subs(curExpr,{N,r,DU}, {5,0.5,DUVec})));
        figure;plot(curExpr_sim);
        curExpr_DU          = diff(curExpr,DU);
        curExpr_DU_abs2     = simplify(rewrite(expand(curExpr_DU),'sincos'));
        pretty(curExpr_DU_abs2);
    end
end
end