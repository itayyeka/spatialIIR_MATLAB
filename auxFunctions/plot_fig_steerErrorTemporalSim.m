function [] = plot_fig_steerErrorTemporalSim()
simCfg              = [];
simCfg.nSensors     = 3;
simCfg.nIterations  = 100;
simCfg.inputFreq    = 10e9;

defaultSimOut   = spatialIIR_singleFreq();
targetAngleVec  = defaultSimOut.targetAngleVec;
N               = defaultSimOut.cfg.nSensors;
duVec           = pi*(cos(targetAngleVec)-cos(defaultSimOut.cfg.thetaS));
nTheta          = defaultSimOut.cfg.nTheta;
if true
    %% subfig_steerErrorTemporalSim_r0
    simCfg.r                = 0;
    theoryBp_dbAbs_norm     = f_theoryBp_dbAbs_norm(duVec,simCfg.r,N);
    simOut                  = spatialIIR_singleFreq(simCfg);
    simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);    
    figure;
    hold on;
    plot(targetAngleVec,theoryBp_dbAbs_norm(:),'*r-','MarkerIndices',1:4:nTheta);
    plot(targetAngleVec,simBp_dbAbs_norm(:),'squareb:','MarkerIndices',3:4:nTheta);
    %% subfig_steerErrorTemporalSim_r06
    simCfg.r                = 0.6;
    theoryBp_dbAbs_norm     = f_theoryBp_dbAbs_norm(duVec,simCfg.r,N);
    simOut                  = spatialIIR_singleFreq(simCfg);
    simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);    
    figure;
    hold on;
    plot(targetAngleVec,theoryBp_dbAbs_norm(:),'*r-','MarkerIndices',1:4:nTheta);
    plot(targetAngleVec,simBp_dbAbs_norm(:),'squareb:','MarkerIndices',3:4:nTheta);
    %% subfig_steerErrorTemporalSim_r08
    simCfg.r                = 0.8;
    theoryBp_dbAbs_norm     = f_theoryBp_dbAbs_norm(duVec,simCfg.r,N);
    simOut                  = spatialIIR_singleFreq(simCfg);
    simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);    
    figure;
    hold on;
    plot(targetAngleVec,theoryBp_dbAbs_norm(:),'*r-','MarkerIndices',1:4:nTheta);
    plot(targetAngleVec,simBp_dbAbs_norm(:),'squareb:','MarkerIndices',3:4:nTheta);
    %% subfig_steerErrorTemporalSim_r095
    simCfg.r                = 0.95;
    theoryBp_dbAbs_norm     = f_theoryBp_dbAbs_norm(duVec,simCfg.r,N);
    simOut                  = spatialIIR_singleFreq(simCfg);
    simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);    
    figure;
    hold on;
    plot(targetAngleVec,theoryBp_dbAbs_norm(:),'*r-','MarkerIndices',1:4:nTheta);
    plot(targetAngleVec,simBp_dbAbs_norm(:),'squareb:','MarkerIndices',3:4:nTheta);
end
end