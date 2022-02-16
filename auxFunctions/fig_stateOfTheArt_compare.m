function [] = fig_stateOfTheArt_compare(simCfg)
    try
        simCfg;
        try
            rCompare   = simCfg.rCompare;
        catch
        end
        try
            rangeError = simCfg.rangeError;
        catch
        end
    catch
        close all;
        simCfg              = f_get_dualFreqCfg();
        rCompare            = 0.5;
        rangeError          = 0;%m
        simCfg.nSensors     = 10;
        simCfg.snr          = inf;
        simCfg.noiseType    = "thermal";
        simCfg.en_mvdr      = 0;
        simCfg.en_CB        = 1;
        simCfg.en_SF_ideal  = 0;
        simCfg.en_SF        = 0;
        simCfg.en_DF        = 1;
    end

    simCfg_base                         = simCfg;
    defaultSimOut                       = spatialIIR_singleFreq();
    targetAngleVec                      = -90 + defaultSimOut.cfg.targetAngleVec*180/pi;
    N                                   = defaultSimOut.cfg.nSensors;
    duVec                               = f_gen_duVec(defaultSimOut.cfg);
    nTheta                              = defaultSimOut.cfg.nTheta;
    simCfg_base.rangeErrorToLambdaRatio = rangeError/defaultSimOut.cfg.lambda;
    if true
        enableMvdr = 0;
        %% MVDR
        if simCfg.en_mvdr
            fs = 8000;
            t = (0:1/fs:1).';
            x1 = cos(2*pi*t*0);
            array = phased.ULA('NumElements', simCfg.nSensors, 'ElementSpacing', simCfg.D);
            array.Element.FrequencyRange = [0.6 1.4]*simCfg.sigFreq;
            fc = simCfg.sigFreq;
            x = collectPlaneWave(array,[x1],[0 0]',fc);
            noise = (10^(-simCfg.snr/10))*((randn(size(x)) + 1i*randn(size(x))))/sqrt(2);
            estimator = phased.MVDREstimator('SensorArray',array,...
                'OperatingFrequency',fc,'DOAOutputPort',true,'NumSignals',1);
            [y,doas] = estimator(x + noise);
            MVDR_lh = plotSpectrum(estimator);
        end
        figure;
        hold on;
        %% CB
        if simCfg.en_CB
            simCfg                  = simCfg_base;
            simCfg.r                = 0;
            simOut                  = spatialIIR_singleFreq(simCfg);
            simBp_dbAbs_norm        = calc_simBp_dbAbs_norm(simOut);
            plot(targetAngleVec,simBp_dbAbs_norm(:));
        end
        %% sf_ideal
        if simCfg.en_SF_ideal
            simCfg                          = simCfg_base;
            simCfg.r                        = rCompare;
            simCfg.snr                      = inf;
            simCfg.rangeErrorToLambdaRatio  = 0;
            simOut                          = spatialIIR_singleFreq(simCfg);
            simBp_dbAbs_norm                = calc_simBp_dbAbs_norm(simOut);
            plot(targetAngleVec,simBp_dbAbs_norm(:));
        end
        %% sf
        if simCfg.en_SF
            simCfg              = simCfg_base;
            simCfg.r            = rCompare;
            simOut              = spatialIIR_singleFreq(simCfg);
            simBp_dbAbs_norm    = calc_simBp_dbAbs_norm(simOut);
            plot(targetAngleVec,simBp_dbAbs_norm(:));
        end
        %% df_rCompare
        if simCfg.en_DF
            simCfg                  = simCfg_base;
            simCfg.r                = 0.5;
            simCfg.kappa            = rCompare;
            simCfg.dF_ratio         = 0.01;
            simOut                  = spatialIIR_dualFreq(simCfg);
            plot(targetAngleVec,simOut.hTwoFreq_err0_dbAbs2_norm(:));
        end
        %% MVDR
        if simCfg.en_mvdr
            MVDR_az = MVDR_lh.XData;
            MVDR_bp = MVDR_lh.YData;
            plot(MVDR_az,MVDR_bp);
        end
    end
    %     fixfig();
end