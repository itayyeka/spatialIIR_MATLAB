function [] = crlbSimpleSim_processSimOut(cfgIn)
close all;
%% path handling
[funcpath, ~, ~]    = fileparts(mfilename('fullpath'));

%% configure
defautCfg.dummy     = 1;

cfgFields = fieldnames(defautCfg);

for cfgFieldID = 1 : numel(cfgFields)
    curCfgField     = cfgFields{cfgFieldID};
    try
        cmdString   = [curCfgField '=cfgIn.(''' curCfgField ''');'];
        eval(cmdString);
    catch
        cmdString   = [curCfgField '=defautCfg.(''' curCfgField ''');'];
        eval(cmdString);
    end
    finalCfg.(curCfgField) = eval([curCfgField ';']);
end

%% process
if true
    %% fetch no-feedback reference results
    simOutDir   = fullfile(funcpath, '..', 'simOut', 'crlbSimpleSim');
    filename    = 'crlbSimpleSimOut_nSensors_3_c_300000kmSec_fSample_20000Mhz_simplifiedSif_noFeedbackREF';
    simOut      = load(fullfile(simOutDir,filename));
    
    try
        bandwidthVec        = simOut.crlbSimpleSim_output.cfg.prmVec.bandwidthVec;
        backoffFactorVec    = simOut.crlbSimpleSim_output.cfg.prmVec.backoffFactorVec;
        rErrFullVec         = simOut.crlbSimpleSim_output.cfg.prmVec.rErrVec;
        thetaFullVec        = simOut.crlbSimpleSim_output.cfg.prmVec.thetaVec;
    catch
        bandwidthVec        = unique(cellfun(@(CELL) CELL.paramSet.bandwidth,       simOut.crlbSimpleSim_output.simOutData_CELL));
        backoffFactorVec    = unique(cellfun(@(CELL) CELL.paramSet.backoffFactor,   simOut.crlbSimpleSim_output.simOutData_CELL));
        rErrFullVec         = unique(cellfun(@(CELL) CELL.paramSet.rErr,            simOut.crlbSimpleSim_output.simOutData_CELL));
        thetaFullVec        = unique(cellfun(@(CELL) CELL.paramSet.theta,           simOut.crlbSimpleSim_output.simOutData_CELL));
    end
    
    if true
        %% CRLB(theta) vs rErr bandwidth = max (REF)
        crlbSimpleSim_fetchSimResults_cfgIn.constantParams.names    = {'backoffFactor',       'theta',          };
        crlbSimpleSim_fetchSimResults_cfgIn.constantParams.values   = {backoffFactorVec(1),   thetaFullVec(1)   };
        crlbSimpleSim_fetchSimResults_cfgIn.simOut                  = simOut;
        [fetchedSimOut] = crlbSimpleSim_fetchSimResults(crlbSimpleSim_fetchSimResults_cfgIn);
        if true
            %% fetch parameters vec
            rErrVec             = cellfun(@(CELL) CELL.paramSet.rErr,       fetchedSimOut);
            crlbThetaVec_REF    = cellfun(@(CELL) CELL.Jinv(1,1),           fetchedSimOut);
        end
    end
    %% fetch simOut
    filename    = 'crlbSimpleSimOut_nSensors_3_c_300000kmSec_fSample_20000Mhz_simplifiedSif.mat';
    simOut      = load(fullfile(simOutDir,filename));
    
    try
        bandwidthVec        = simOut.crlbSimpleSim_output.cfg.prmVec.bandwidthVec;
        backoffFactorVec    = simOut.crlbSimpleSim_output.cfg.prmVec.backoffFactorVec;
        rErrFullVec         = simOut.crlbSimpleSim_output.cfg.prmVec.rErrVec;
        thetaFullVec        = simOut.crlbSimpleSim_output.cfg.prmVec.thetaVec;
    catch
        bandwidthVec        = unique(cellfun(@(CELL) CELL.paramSet.bandwidth,       simOut.crlbSimpleSim_output.simOutData_CELL));
        backoffFactorVec    = unique(cellfun(@(CELL) CELL.paramSet.backoffFactor,   simOut.crlbSimpleSim_output.simOutData_CELL));
        rErrFullVec         = unique(cellfun(@(CELL) CELL.paramSet.rErr,            simOut.crlbSimpleSim_output.simOutData_CELL));
        thetaFullVec        = unique(cellfun(@(CELL) CELL.paramSet.theta,           simOut.crlbSimpleSim_output.simOutData_CELL));
    end
    
    
    %{
    bandwidth
    backoffFactor
    rErr
    theta
    %}
    %% CRLB(theta) vs bandwidth rErr = 0
    crlbSimpleSim_fetchSimResults_cfgIn                         = [];
    refBackoffFactor                                            = backoffFactorVec(1);
    crlbSimpleSim_fetchSimResults_cfgIn.constantParams.names    = {'backoffFactor',     'theta',            'rErr'  };
    crlbSimpleSim_fetchSimResults_cfgIn.constantParams.values   = {refBackoffFactor,    thetaFullVec(1),    0       };
    crlbSimpleSim_fetchSimResults_cfgIn.simOut                  = simOut;
    [fetchedSimOut] = crlbSimpleSim_fetchSimResults(crlbSimpleSim_fetchSimResults_cfgIn);
    if true
        %% fetch parameters vec
        bandwidthVec    = cellfun(@(CELL) CELL.paramSet.bandwidth,  fetchedSimOut);
        crlbThetaVec    = cellfun(@(CELL) CELL.Jinv(1,1),           fetchedSimOut);
        %% generate MESH
        figure;
        plot(bandwidthVec,db(abs([reshape(crlbThetaVec_REF(1:length(crlbThetaVec)),[],1) crlbThetaVec(:)])));
        xlabel('signal bandwidth [normalized]');
        ylabel('CRLB_\theta [dB]');
        legend({'no-feedback',['with feedback backoff = ' num2str(refBackoffFactor) ', rErr = ' num2str(0)]})
    end
    
    %% CRLB(theta) vs (rErr & bandwidth)
    plotBackoffFactor                                           = backoffFactorVec(2);
    crlbSimpleSim_fetchSimResults_cfgIn                         = [];
    crlbSimpleSim_fetchSimResults_cfgIn.constantParams.names    = {'backoffFactor',         'theta'         };
    crlbSimpleSim_fetchSimResults_cfgIn.constantParams.values   = {plotBackoffFactor,   thetaFullVec(1) };
    crlbSimpleSim_fetchSimResults_cfgIn.simOut                  = simOut;
    [fetchedSimOut] = crlbSimpleSim_fetchSimResults(crlbSimpleSim_fetchSimResults_cfgIn);
    if true
        %% fetch parameters vec
        rErrVec         = cellfun(@(CELL) CELL.paramSet.rErr,       fetchedSimOut);
        bandwidthVec    = cellfun(@(CELL) CELL.paramSet.bandwidth,  fetchedSimOut);
        crlbThetaVec    = cellfun(@(CELL) CELL.Jinv(1,1),           fetchedSimOut);
        %% generate MESH
        figure;
        trisurf(delaunay(rErrVec,bandwidthVec),rErrVec,bandwidthVec,db(abs(crlbThetaVec)));
        xlabel('r error [m]');
        ylabel('signal bandwidth [normalized]');
        zlabel('CRLB_\theta');
        title({...
            ['CRLB_{\theta} vs. rErr and signal bandwidth']...
            ['rErr up to 2*\lambda (\lambda =' num2str(rErrVec(end)/2) ')']...
            ['backoff factor = ' num2str(plotBackoffFactor)] ...
            })
    end
    %% CRLB(range) vs (rErr & bandwidth)
    if false
        %% fetch parameters vec
        rErrVec         = cellfun(@(CELL) CELL.paramSet.rErr,       fetchedSimOut);
        bandwidthVec    = cellfun(@(CELL) CELL.paramSet.bandwidth,  fetchedSimOut);
        crlbThetaVec    = cellfun(@(CELL) CELL.Jinv(2,2),           fetchedSimOut);
        %% generate MESH
        figure;
        trisurf(delaunay(rErrVec,bandwidthVec),rErrVec,bandwidthVec,db(abs(crlbThetaVec)));
        xlabel('r error [m]');
        ylabel('signal bandwidth [normalized]');
        zlabel('CRLB_{range}');
    end
end
end