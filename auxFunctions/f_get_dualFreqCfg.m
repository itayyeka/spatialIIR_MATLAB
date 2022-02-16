function [simCfg] = f_get_dualFreqCfg(simCfg)
    try
        simCfg;
    catch
        simCfg = [];
    end
    simCfg = f_get_singleFreqCfg(simCfg);    
    try
        dF_ratio  = simCfg.dF_ratio;
    catch
        dF_ratio  = 1e-1;
    end
    simCfg.dF_ratio = dF_ratio;    
    try
        kappa  = simCfg.kappa;
    catch
        kappa  = simCfg.r;
    end
    simCfg.kappa = kappa;
end