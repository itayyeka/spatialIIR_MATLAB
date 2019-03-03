function [] = spatialIIR_symbolicTheroyWrapper()
clear;
clc;
close all;

%% path handling
[funcPath, ~, ~] = fileparts(mfilename('fullpath'));
try
    spatialIIR_MATLAB_subFunctions_indicator;
catch
    addpath(genpath(fullfile(funcPath,'subFunctions')));
    spatialIIR_MATLAB_subFunctions_indicator;
end

%% configure
scriptFlags.crlbSimpleSim = 0;
scriptFlags.beamwidthCalc = 1;
scriptFlags.tauErrorEffect = 0;

%% execute
beamwidthMAT = [];
for nSensors = 3 : 4
    spatialIIR_symbolicTheroy_cfgIn.nSensors    = nSensors;
    
    simOut  = spatialIIR_symbolicTheroy(spatialIIR_symbolicTheroy_cfgIn,scriptFlags);
    
    simOut.beamwidthCalc_output.beamwidth_withFeedback;
    
%     beamwidthMAT = [beamwidth_withFeedback simOut.beamwidthCalc_output.beamwidth_withFeedback];
%     
%     rPoleVec
%     legendSTR = cellfun(@(rPole) ['nSensors = '], num2str());
end
end