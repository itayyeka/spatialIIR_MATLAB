function [] = plot_fig_stabilization()
    simCfg              = [];
    simCfg.nSensors     = 3;
    simCfg.nIterations  = 100;
    simCfg.inputFreq    = 10e9;
    simCfg.nTheta       = 1;
    
    simCfg_base     = simCfg;
    defaultSimOut   = spatialIIR_singleFreq(simCfg);
    targetAngleVec  = defaultSimOut.targetAngleVec;
    N               = defaultSimOut.cfg.nSensors;
    duVec           = pi*(cos(targetAngleVec)-cos(defaultSimOut.cfg.thetaS));
    nTheta          = simCfg.nTheta;
    markerIndices   = 1:15:nTheta;
    %% stability_vs_r
    if false
        %% subfig_steerErrorTemporalSim_r06
        simCfg                  = simCfg_base;
        simCfg.r                = 0.6;
        simOut                  = spatialIIR_singleFreq(simCfg);
        fId = figure;
        iterId = (1 : length(simOut.hMat(:))) / simOut.cfg.nSamplesIter;
        plot(iterId, db(abs(simOut.hMat(:))))
        hold on
        %% subfig_steerErrorTemporalSim_r09
        simCfg                  = simCfg_base;
        simCfg.r                = 0.9;
        simOut                  = spatialIIR_singleFreq(simCfg);
        plot(iterId, db(abs(simOut.hMat(:))))
        %% subfig_steerErrorTemporalSim_r1
        simCfg                  = simCfg_base;
        simCfg.r                = 1;
        simOut                  = spatialIIR_singleFreq(simCfg);
        plot(iterId, db(abs(simOut.hMat(:))))
        %% subfig_steerErrorTemporalSim_r11
        simCfg                  = simCfg_base;
        simCfg.r                = 1.1;
        simOut                  = spatialIIR_singleFreq(simCfg);
        plot(iterId, db(abs(simOut.hMat(:))))
        %% cfg
        ylim([0 50]);
    end
    %% stability duration - r
    if false
        %% cfg
        nr = 100;
        %% sim
        rVec = linspace(0,.99,nr);
        stblDur = zeros(size(rVec));
        finalVal = zeros(size(rVec));
        for r = rVec
            %% temporal sim
            simCfg              = simCfg_base;
            simCfg.r            = r;
            simCfg.nIterations  = 1000;
            simOut              = spatialIIR_singleFreq(simCfg);
            %% process
            finalAbs = abs(simOut.hMat(end));
            finalVal(rVec == r) = finalAbs;
            try
                stblDur(rVec == r) = find(abs(simOut.hMat)>0.99*finalAbs, true, 'first');
            catch
            end
        end
        %% plot
        fId = figure;
        subplot(1,4,1:2);
        stblIterId = stblDur / simOut.cfg.nSamplesIter;
        plot(rVec, stblIterId);
        ylim([0 100]);        
        subplot(1,4,3:4);
        plot(rVec, db(finalVal));
        fixfig();
    end
    %% stability duration - N
    if true
        %% cfg
        minN = 2;
        maxN = 10;
        %% sim        
        nVec = minN:maxN;
        stblDur = zeros(size(nVec));
        finalVal = zeros(size(nVec));
        for N = nVec
            nId = N-minN+1;
            %% temporal sim
            simCfg              = simCfg_base;
            simCfg.r            = 0.6;
            simCfg.nSensors     = N;
            simCfg.nIterations  = 1000;
            simOut              = spatialIIR_singleFreq(simCfg);
            %% process
            finalAbs = abs(simOut.hMat(end));
            finalVal(nId) = finalAbs;
            try
                stblDur(nId) = find(abs(simOut.hMat)>0.99*finalAbs, true, 'first');
            catch
            end
        end
        %% plot
        fId = figure;
        subplot(1,4,1:2);
        stblIterId = stblDur / simOut.cfg.nSamplesIter;
        plot(nVec, stblIterId);
        ylim([0 100]);        
        subplot(1,4,3:4);
        plot(nVec, db(finalVal));
        fixfig();
    end
end