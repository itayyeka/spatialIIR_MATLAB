function [] = plot_fig_rangError()
simCfg              = [];
simCfg.nSensors     = 3;
simCfg.nIterations  = 100;
simCfg.inputFreq    = 10e9;
simCfg.nTheta       = 300;

simCfg_base     = simCfg;
defaultSimOut   = spatialIIR_singleFreq(simCfg);
targetAngleVec  = defaultSimOut.targetAngleVec;
N               = defaultSimOut.cfg.nSensors;
duVec           = pi*(cos(targetAngleVec)-cos(defaultSimOut.cfg.thetaS));
nTheta          = simCfg.nTheta;
markerIndices   = 1:15:nTheta;
if true    
    %% subfig_steerErrorTemporalSim_r06_err0
    if true
        simCfg                  = simCfg_base;
        simCfg.r                = 0.6;
        theoryBp_dbAbs_norm     = f_theoryBp_dbAbs_norm(duVec,simCfg.r,N);
        simOut                  = spatialIIR_singleFreq(simCfg);
        simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);
        figure;
        hold on;
        plot(duVec/pi,theoryBp_dbAbs_norm(:),'*b-','MarkerIndices',markerIndices);
        %plot(duVec,simBp_dbAbs_norm(:),'og:','MarkerIndices',3:8:nTheta);
    end
    %% subfig_steerErrorTemporalSim_r06_err01
    if true
        simCfg                          = simCfg_base;
        simCfg.r                        = 0.6;
        simCfg.rangeErrorToLambdaRatio  = 0.1;
        simOut                          = spatialIIR_singleFreq(simCfg);
        simBp_dbAbs_norm                = calc_simBp_dbAbs_norm(simOut);        
        plot(duVec/pi,simBp_dbAbs_norm(:),'squarer--','MarkerIndices',markerIndices);
    end
    %% subfig_steerErrorTemporalSim_r06_err03
    if true
        simCfg                          = simCfg_base;
        simCfg.r                        = 0.6;
        simCfg.rangeErrorToLambdaRatio  = 0.3;
        simOut                          = spatialIIR_singleFreq(simCfg);
        simBp_dbAbs_norm                = calc_simBp_dbAbs_norm(simOut);        
        plot(duVec/pi,simBp_dbAbs_norm(:),'diamondr-','MarkerIndices',markerIndices);
    end
    %% cfg
    ylim([-50 0]);
end
end