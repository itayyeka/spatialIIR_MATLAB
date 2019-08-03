function [] = plot_fig_dualfreq_rangeErrorHighSnr()
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
    simCfg.rangeErrorToLambdaRatio  = 0.1;
    simCfg.thetaS                   = pi/2;
    simCfg.snr                      = inf;
    simCfg.kappa                    = simCfg.r;
    [simOut_range_rangeErr01]       = spatialIIR_dualFreq(simCfg);
    %% range error
    simCfg.rangeErrorToLambdaRatio  = 0.3;
    [simOut_range_rangeErr03]       = spatialIIR_dualFreq(simCfg);
end
%% plot
nTheta          = simOut_range_rangeErr01.nTheta;
targetAngleVec  = simOut_range_rangeErr03.targetAngleVec/pi;
indicesDistance = 2;
indicesVec1     = 1:indicesDistance*4:nTheta;
indicesVec2     = 1+indicesDistance:indicesDistance*4:nTheta;
indicesVec3     = 1+2*indicesDistance:indicesDistance*4:nTheta;
figure;
subplot(1,2,1);
hold on;
plot(targetAngleVec,simOut_range_rangeErr01.theoryBp_dbAbs_norm(:),'*b-','MarkerIndices',indicesVec1);
plot(targetAngleVec,simOut_range_rangeErr01.bp_singleReq_norm,'squarer-','MarkerIndices',indicesVec2);
plot(targetAngleVec,simOut_range_rangeErr01.hTwoFreq_err0_dbAbs2_norm(:),'diamondg-','MarkerIndices',indicesVec3);
subplot(1,2,2);
hold on;
plot(targetAngleVec,simOut_range_rangeErr03.theoryBp_dbAbs_norm(:),'*b-','MarkerIndices',indicesVec1);
plot(targetAngleVec,simOut_range_rangeErr03.bp_singleReq_norm,'squarer-','MarkerIndices',indicesVec2);
plot(targetAngleVec,simOut_range_rangeErr03.hTwoFreq_err0_dbAbs2_norm(:),'diamondg-','MarkerIndices',indicesVec3);
end