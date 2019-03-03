clear all;
close all;
clc;

try
    auxFunctionsIndicatorFunctions;
catch
    restoredefaultpath;
    [funcPath,~,~]  = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(funcPath,'auxFunctions')));
    auxFunctionsIndicatorFunctions;
end

%% general case BP (\label{eqn_generalCaseBp})
if true
    %% symbolics
    syms N DU r dPhi positive;
    
    %% generating basic terms
    ad  = (1/N)*(1-exp(-1i*N*DU))/(1-exp(-1i*DU));
    bd  = ad;
    
    %% full expression transfer function
    hDen    = simplify(1-r*ad*exp(1i*dPhi));
    h       = simplify(bd/hDen);
    hAbs2   = simplify(rewrite(expand(h.*conj(h)),'sincos'));
    
    %% hOpt
    hAbs2_dPhi0         = subs(hAbs2, dPhi, 0);
    [lim_hAbs2,report]  = get_funcLimit(hAbs2_dPhi0,DU,0);
    hAbs2Rel            = simplify(rewrite(expand(hAbs2/lim_hAbs2),'sincos'));
    pretty(hAbs2Rel)
end

%% Figures : fig_feedbackULA_HPBW_Nx_vs_N_variousR,fig_feedbackULA_beamwidth_limit_r_dependent
if false
    plot_fig_HPBW(hAbs2Rel);
end

%% Figure : fig_firstSidelobeGain_CB
plot_fig_sideLobes(hAbs2Rel);
