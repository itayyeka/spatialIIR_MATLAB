function [duVec] = f_gen_duVec(simCfg)
    try
        duMin   = f_convert_angle2du(simCfg.thetaMin,simCfg);
        duMax   = f_convert_angle2du(simCfg.thetaMax,simCfg);
    catch
        try
            duMin   = simCfg.duMin;
            duMax   = simCfg.duMax;
        catch
            duMin   = f_convert_angle2du(0,simCfg);
            duMax   = f_convert_angle2du(pi,simCfg);
        end
    end
    duLimits    = sort([duMin,duMax]);
    duVec       = linspace(duLimits(1),duLimits(2),simCfg.nTheta);
end

