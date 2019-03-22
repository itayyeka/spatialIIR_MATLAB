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
DUVec       = linspace(-0.05,0.05,200);
hSteerErr   = subs(h,dPhi,0);
if true
    %% sidelobes
    if false
        hNum            = (1-exp(1i*N*DU))/(N*(1-exp(1i*DU)));
        hDen            = 1-r*hNum;
        hNumAbs2        = simplify(rewrite(expand(hNum*conj(hNum)),'sincos'));
        hDenAbs2        = simplify(rewrite(expand(hDen*conj(hDen)),'sincos'));
        steerH          = hNum/hDen;
        steerHAbs2      = simplify(rewrite(expand(steerH*conj(steerH)),'sincos'));
        pretty(hDenAbs2);
        DUVec           = linspace(-pi,pi,1000);
        NVal            = 10;
        rVal            = 0.9;
        hNumVal         = subs(hNumAbs2,{DU,N,r},{DUVec, NVal, rVal});
        hDenVal         = subs(hDenAbs2,{DU,N,r},{DUVec, NVal, rVal});
        steerHAbs2Val   = subs(steerHAbs2,{DU,N,r},{DUVec, NVal, rVal});
        
        hNumVal_norm        = hNumVal/max(abs(hNumVal));
        hDenVal_norm        = hDenVal/max(abs(hDenVal));
        steerHAbs2Val_norm  = steerHAbs2Val/max(abs(steerHAbs2Val));
        
        figure;plot(DUVec,db([hNumVal_norm(:) hDenVal_norm(:) steerHAbs2Val_norm(:)]));
        hDiff           = simplify(rewrite(expand(diff(hDenAbs2,DU)),'exp'));
        pretty(hDiff);
        
    end
    %% HPBW
    syms x P;
    nTaylor     = 4;
    taylorPoly  = 1;
    nVal        = 1000;
    rVal        = 0.2;
    PVal        = 1/2;
    if true
        D                   = sin(x)/x;
        fullExpr            = (P*r^2-(1-r)^2)*D^2-2*P*r*cos(x)*D+P;
        fullExpr_taylor     = get_funcLimit(fullExpr,x,0);
        curExpr             = fullExpr;
        for orderId = 1 : nTaylor
            curExpr             = diff(curExpr,x);
            curDiffLim          = get_funcLimit(curExpr,x,0);
            fullExpr_taylor     = fullExpr_taylor + curDiffLim*x^orderId/factorial(orderId);
        end
        xVec                                = linspace(1e-6,3,100);
        fullExpr_val                        = eval(subs(fullExpr,{P r N x},{PVal rVal nVal xVec}));
        fullExpr_taylor_val                 = eval(subs(fullExpr_taylor,{P r N x},{PVal rVal nVal xVec}));
        figure;plot(xVec,[fullExpr_val(:) fullExpr_taylor_val(:)]);
        fullExprTaylorCoeffs                = simplify(coeffs(fullExpr_taylor,x,'All'));
        fullExprTaylorCoeffs_HPBW           = simplify(subs(fullExprTaylorCoeffs,P,PVal));
        fullExprTaylorCoeffs_HPBW_roots     = roots(fullExprTaylorCoeffs_HPBW);
        fullExprTaylorCoeffs_HPBW_rootsVec  = eval(subs(fullExprTaylorCoeffs_HPBW_roots,{P,r},{1/2 0.5}));
        rootId = ...
            find(...
            cellfun(@isreal,num2cell(fullExprTaylorCoeffs_HPBW_rootsVec)) ...
            .* ...
            cellfun(@(x) real(x)>0 ,num2cell(fullExprTaylorCoeffs_HPBW_rootsVec)) ...
            );
        fullExprTaylorCoeffs_HPBW_root          = simplify(expand(fullExprTaylorCoeffs_HPBW_roots(rootId)));
        rVec                                    = linspace(0,1,100);
        fullExprTaylorCoeffs_HPBW_rootValVec    = eval(subs(fullExprTaylorCoeffs_HPBW_root,r,rVec));
        figure;plot(rVec,fullExprTaylorCoeffs_HPBW_rootValVec);
    end
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