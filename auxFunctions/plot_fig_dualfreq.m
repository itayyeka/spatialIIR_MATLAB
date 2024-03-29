function [] = plot_fig_dualfreq()
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
    simCfg.snr                      = 10;
    simCfg.kappa                    = simCfg.r;
    [simOut_range_snrHigh]          = spatialIIR_dualFreq(simCfg);
    %% range error
    simCfg.snr                      = 3;
    [simOut_range_snrLow]           = spatialIIR_dualFreq(simCfg);
end
%% plot
nTheta          = simOut_range_snrHigh.cfg.nTheta;
targetAngleVec  = simOut_range_snrHigh.cfg.targetAngleVec/pi;
indicesDistance = 2;
indicesVec1     = 1:indicesDistance*4:nTheta;
indicesVec2     = 1+indicesDistance:indicesDistance*4:nTheta;
indicesVec3     = 1+2*indicesDistance:indicesDistance*4:nTheta;
figure;
subplot(1,2,1);
hold on;
plot(targetAngleVec,simOut_range_snrHigh.theoryBp_dbAbs_norm(:),'*b-','MarkerIndices',indicesVec1);
plot(targetAngleVec,simOut_range_snrHigh.bp_singleReq_norm,'squarer-','MarkerIndices',indicesVec2);
plot(targetAngleVec,simOut_range_snrHigh.hTwoFreq_err0_dbAbs2_norm(:),'diamondg-','MarkerIndices',indicesVec3);
subplot(1,2,2);
hold on;
plot(targetAngleVec,simOut_range_snrLow.theoryBp_dbAbs_norm(:),'*b-','MarkerIndices',indicesVec1);
plot(targetAngleVec,simOut_range_snrLow.bp_singleReq_norm,'squarer-','MarkerIndices',indicesVec2);
plot(targetAngleVec,simOut_range_snrLow.hTwoFreq_err0_dbAbs2_norm(:),'diamondg-','MarkerIndices',indicesVec3);
end