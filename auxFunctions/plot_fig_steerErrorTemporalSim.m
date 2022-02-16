function [] = plot_fig_steerErrorTemporalSim(simCfg)
    try
        simCfg;
    catch
        close all;
        simCfg                          = [];
        simCfg.nSensors                 = 3;
        simCfg.nIterations              = 100;
        simCfg.inputFreq                = 10e9;
        simCfg.rangeErrorToLambdaRatio  = 0;
        simCfg.snr                      = 5;%10^(-snr/10)
    end

    simCfg_base     = simCfg;
    defaultSimOut   = spatialIIR_singleFreq();
    targetAngleVec  = defaultSimOut.cfg.targetAngleVec;
    N               = defaultSimOut.cfg.nSensors;
    duVec           = f_gen_duVec(defaultSimOut.cfg);
    nTheta          = defaultSimOut.cfg.nTheta;
    if true
        %% MVDR
        fs = 8000;
        t = (0:1/fs:1).';
        x1 = cos(2*pi*t*0);
        array = phased.ULA('NumElements',defaultSimOut.cfg.nSensors,'ElementSpacing',defaultSimOut.cfg.D);
        array.Element.FrequencyRange = [0.6 1.4]*simCfg.inputFreq;
        fc = simCfg.inputFreq;
        x = collectPlaneWave(array,[x1],[0 0]',fc);
        snr_amp = 10^(-simCfg.snr/10);
        noise = snr_amp*(randn(size(x)) + 1i*randn(size(x)))/sqrt(2);
        estimator = phased.MVDREstimator('SensorArray',array,...
            'OperatingFrequency',fc,'DOAOutputPort',true,'NumSignals',1);
        [y,doas] = estimator(x + noise);
        doas = broadside2az(sort(doas),[20 -5]);        
        specMVDR_lh = plotSpectrum(estimator);
        figure;
        specMVDR_az = specMVDR_lh.XData;
        specMVDR_dB = specMVDR_lh.YData;
        plot(specMVDR_az,specMVDR_dB);
        %% bartlett
        if true
            simCfg                  = simCfg_base;
            simCfg.r                = 0;
            theoryBp_dbAbs_norm     = f_theoryBp_dbAbs_norm(duVec,simCfg);
            simOut                  = spatialIIR_singleFreq(simCfg);
            simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);            
            hold on;
            plot(-90+targetAngleVec*180/pi,simBp_dbAbs_norm(:));
        end
        %% r06
        if true
            simCfg                  = simCfg_base;
            simCfg.r                = 0.6;
            theoryBp_dbAbs_norm     = f_theoryBp_dbAbs_norm(duVec,simCfg);
            simOut                  = spatialIIR_singleFreq(simCfg);
            simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);            
            hold on;
            plot(-90+targetAngleVec*180/pi,simBp_dbAbs_norm(:));
        end
        %% DF_k06
        if true
            simCfg                  = simCfg_base;
            simCfg.r                = 0.6;
            theoryBp_dbAbs_norm     = f_theoryBp_dbAbs_norm(duVec,simCfg);
            simOut                  = spatialIIR_singleFreq(simCfg);
            simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);            
            hold on;
            plot(-90+targetAngleVec*180/pi,simBp_dbAbs_norm(:));
        end
        %% subfig_steerErrorTemporalSim_r09
        if true
            simCfg                  = simCfg_base;
            simCfg.r                = 0.9;
            theoryBp_dbAbs_norm     = f_theoryBp_dbAbs_norm(duVec,simCfg);
            simOut                  = spatialIIR_singleFreq(simCfg);
            simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);            
            hold on;
            plot(-90+targetAngleVec*180/pi,simBp_dbAbs_norm(:));
        end        
    end
end