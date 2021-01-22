clear all;
close all;
clc;

try
    auxFunctionsIndicatorFunctions;
catch
    %%restoredefaultpath;
    [funcPath,~,~]  = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(funcPath,'auxFunctions')));
    auxFunctionsIndicatorFunctions;
end

%% general case BP (\label{eqn_generalCaseBp})
syms dPhi DU N;
if false
    [hAbs2Rel]          = H_development();
    hAbs2Rel_steerError = subs(hAbs2Rel,dPhi,0);
end
%% Figures : fig_feedbackULA_HPBW_Nx_vs_N_variousR,fig_feedbackULA_beamwidth_limit_r_dependent
if false
    plot_fig_HPBW(hAbs2Rel);
end
%% sidelobes
if false    
    %% symbolics
    dPhi_sideLobes      = 0;
    DU_sideLobes        = 3*pi/N;
    hAbs2Rel_sidelobes  = subs(hAbs2Rel,{dPhi DU},{dPhi_sideLobes,DU_sideLobes});
    %% plot fig_firstSidelobeGain_CB
    plot_fig_sideLobes(hAbs2Rel);
end
%% directivity
if false
    directivityAnalysis(1);
end

%% Figure : fig_steerErrorTemporalSim
if false
    plot_fig_steerErrorTemporalSim;
end

%% Figure : fig_hDUDTContour
if false
    plot_fig_hDUDTContour(hAbs2Rel);
end

%% Figure : fig_rangError
if false
    plot_fig_rangError;
end

%% application
if false
    plot_fig_dualfreq_perfectAlignHighSnr;
    plot_fig_dualfreq_perfectAlignLowSnr;
    plot_fig_dualfreq_rangeErrorHighSnr;
end