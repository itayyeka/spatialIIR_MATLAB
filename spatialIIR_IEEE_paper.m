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
    hAbs2Rel    = H_development();
end
%% Figures : fig_feedbackULA_HPBW_Nx_vs_N_variousR,fig_feedbackULA_beamwidth_limit_r_dependent
if false
    plot_fig_HPBW(hAbs2Rel);
end

%% Figure : fig_firstSidelobeGain_CB
plot_fig_sideLobes(hAbs2Rel);
