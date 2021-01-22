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

%% Figures : fig_feedbackULA_HPBW_Nx_vs_N_variousR,fig_feedbackULA_beamwidth_limit_r_dependent
if false
    plot_fig_stabilization();
end
if false
    plot_fig_beamThining();
end