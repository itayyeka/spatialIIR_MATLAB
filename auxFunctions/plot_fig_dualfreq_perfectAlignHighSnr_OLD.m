function [] = plot_fig_dualfreq_perfectAlignHighSnr()
%% collect simulations
if true
    %% perfect alignmnet
    simCfg                          = [];
    simCfg.nSensors                 = 3;
    simCfg.nIterations              = 100;
    simCfg.inputFreq                = 10e9;
    simCfg.dF_ratio                 = 1e-1;
    simCfg.targetRange_samples      = 32;
    simCfg.r                        = 0.6;
    simCfg.rangeErrorToLambdaRatio  = 0;
    simCfg.thetaS                   = pi/2;
    simCfg.snr                      = inf;
    simCfg.kappa                    = simCfg.r;
    [simOut_ideal]                  = spatialIIR_dualFreq(simCfg);
    %% range error
    simCfg.thetaS                   = 3*pi/4;
    [simOut_rangeErr]               = spatialIIR_dualFreq(simCfg);
end
%% plot
nTheta          = simOut_ideal.nTheta;
targetAngleVec  = simOut_ideal.targetAngleVec/pi;
indicesDistance = 4;
indicesVec1     = 1:indicesDistance*4:nTheta;
indicesVec2     = 1+indicesDistance:indicesDistance*4:nTheta;
indicesVec3     = 1+2*indicesDistance:indicesDistance*4:nTheta;
figure;
subplot(1,2,1);
hold on;
plot(targetAngleVec,simOut_ideal.theoryBp_dbAbs_norm(:),'*b-','MarkerIndices',indicesVec1);
plot(targetAngleVec,simOut_ideal.bp_singleReq_norm,'squarer-','MarkerIndices',indicesVec2);
plot(targetAngleVec,simOut_ideal.hTwoFreq_err0_dbAbs2_norm(:),'diamondg-','MarkerIndices',indicesVec3);
subplot(1,2,2);
hold on;
plot(targetAngleVec,simOut_rangeErr.theoryBp_dbAbs_norm(:),'*b-','MarkerIndices',indicesVec1);
plot(targetAngleVec,simOut_rangeErr.bp_singleReq_norm,'squarer-','MarkerIndices',indicesVec2);
plot(targetAngleVec,simOut_rangeErr.hTwoFreq_err0_dbAbs2_norm(:),'diamondg-','MarkerIndices',indicesVec3);
end