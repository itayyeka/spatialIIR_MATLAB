function [] = compare_CB_to_DF(nSensors,r,snr)
    simCfg = f_get_dualFreqCfg();
    try
        simCfg.rCompare = r;
    catch
        simCfg.rCompare = 0.5;
    end
    simCfg.rangeError = 0;%m
    try
        simCfg.nSensors = nSensors;
    catch
        simCfg.nSensors = 10;
    end
    try
        simCfg.snr = snr;
    catch
        simCfg.snr = inf;
    end
    simCfg.en_mvdr      = 0;
    simCfg.en_CB        = 1;
    simCfg.en_SF_ideal  = 0;
    simCfg.en_SF        = 0;
    simCfg.en_DF        = 1;
    fig_stateOfTheArt_compare(simCfg);
end