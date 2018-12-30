function [] = spatialIIR_symbolicTheroyWrapper(cfgIn,scriptFlags)
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
try
    cfgIn;
catch
    cfgIn = [];
end
try
    scriptFlags.crlbSimpleSim;
catch
    scriptFlags.crlbSimpleSim = 1;
end
try
    nSensors = cfgIn.nSensors;
catch
    nSensors = 3;
end

%% symbolics
alphaV              = sym('alpha',[nSensors 1]);
alphaT              = transpose(alphaV);
alphaH              = conj(alphaT);
betaV               = sym('beta',[nSensors 1]);
betaT               = transpose(betaV);
betaH               = conj(betaT);
syms                omega;
assume(             omega, 'real');
syms                targetRange;
assume(             targetRange, 'real');
syms                c;
assume(             c, 'real');
tau                 = targetRange/c;
syms                D;
assume(             D, 'real');
syms                theta;
assume(             theta, 'real');
syms                sigmaX;
assume(             sigmaX, 'real');
tauTheta            = D*cos(theta)/c;
sensorID_zeroBased  = 0:(nSensors-1);
sensorID_oneBased   = 1:nSensors;
d                   = reshape(exp(-1i*omega*tauTheta*sensorID_zeroBased),[],1);
dT                  = transpose(d);
dH                  = conj(dT);
g                   = betaT*d*exp(-1i*omega*tau);
A                   = -1i*omega*D*sin(theta)*diag(sensorID_zeroBased)/c;
AT                  = transpose(A);
AH                  = conj(AT);
B                   = d*dT*A-A*d*dT;
BT                  = transpose(B);
BH                  = conj(BT);

crlbVariables.nSensors              = nSensors;
crlbVariables.alphaV                = alphaV;
crlbVariables.alphaT                = alphaT;
crlbVariables.alphaH                = alphaH;
crlbVariables.betaV                 = betaV;
crlbVariables.betaT                 = betaT;
crlbVariables.betaH                 = betaH;
crlbVariables.omega                 = omega;
crlbVariables.tau                   = tau;
crlbVariables.targetRange           = targetRange;
crlbVariables.c                     = c;
crlbVariables.D                     = D;
crlbVariables.theta                 = theta;
crlbVariables.sigmaX                = sigmaX;
crlbVariables.tauTheta              = tauTheta;
crlbVariables.sensorID_zeroBased    = sensorID_zeroBased;
crlbVariables.sensorID_oneBased     = sensorID_oneBased;
crlbVariables.d                     = d;
crlbVariables.dT                    = dT;
crlbVariables.dH                    = dH;
crlbVariables.g                     = g;
crlbVariables.A                     = A;
crlbVariables.AT                    = AT;
crlbVariables.AH                    = AH;
crlbVariables.B                     = B;
crlbVariables.BT                    = BT;
crlbVariables.BH                    = BH;


%% sub scripts
if scriptFlags.crlbSimpleSim
    %% configure
    crlbSimpleSimCfg.backoffFactor_min                      = 1;
    crlbSimpleSimCfg.backoffFactor_max                      = 1;
    crlbSimpleSimCfg.backoffFactor_nValues                  = 1;
    crlbSimpleSimCfg.rErr_min                               = -1;
    crlbSimpleSimCfg.rErr_max                               = 1;
    crlbSimpleSimCfg.rErr_nValues                           = 11;
    crlbSimpleSimCfg.thetaSim_azimuthalWidth                = 0;
    crlbSimpleSimCfg.thetaSim_nPoints                       = 1;
    crlbSimpleSimCfg.polesThetaVec                          = repmat(pi/3,1,nSensors-1);
    crlbSimpleSimCfg.RangeVal                               = 1000;
    crlbSimpleSimCfg.DVal                                   = 0.01;
    crlbSimpleSimCfg.cVal                                   = 3e8;
    crlbSimpleSimCfg.fSample                                = 20e9; 
    crlbSimpleSimCfg.syncSigCfg.duration_SAMPLES            = 1024;
    crlbSimpleSimCfg.syncSigCfg.baseFreq_relative           = 0.2;
    crlbSimpleSimCfg.syncSigCfg.bandwidth_relative_max      = 0.3;
    crlbSimpleSimCfg.syncSigCfg.bandwidth_relative_nValues  = 30;
    crlbSimpleSimCfg.integralNPoints                        = 128;
    %% execute
    crlbSimpleSim_output    = crlbSimpleSim(crlbVariables,crlbSimpleSimCfg);
    finaCfg                 = crlbSimpleSim_output.cfg;
    simOut_savedName        = ['crlbSimpleSimOut_nSensors_' num2str(nSensors) '_c_' num2str(finaCfg.cVal/1000) 'kmSec_fSample_' num2str(finaCfg.fSample/1e6) 'Mhz'];
    savedFilename           = [regexprep(simOut_savedName,'\.','_') '.mat'];
    outputDir               = fullfile(funcPath,'simOut','crlbSimpleSim');
    if ~isdir(outputDir)
        mkdir(outputDir);
    end
    savedFilePath           = fullfile(outputDir,savedFilename);
    save(savedFilePath,'crlbSimpleSim_output');
end
end