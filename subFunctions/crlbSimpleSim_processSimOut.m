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
    %% fetch simOut
    crlbSimpleSim_fetchSimResults_cfgIn.simOutDir               = fullfile(funcpath, '..', 'simOut', 'crlbSimpleSim');
    crlbSimpleSim_fetchSimResults_cfgIn.filename                = 'crlbSimpleSimOut_nSensors_3_c_300000kmSec_fSample_20000Mhz.mat';
    crlbSimpleSim_fetchSimResults_cfgIn.constantParams.names    = {'backoffFactor', 'theta' };
    crlbSimpleSim_fetchSimResults_cfgIn.constantParams.values   = {1,               pi/3    };
    [fetchedSimOut] = crlbSimpleSim_fetchSimResults(crlbSimpleSim_fetchSimResults_cfgIn);
    %% plot #1 CRLB vs (rErr & bandwidth)
    if true
        %% fetch parameters vec
        rErrVec         = cellfun(@(CELL) CELL.paramSet.rErr,       fetchedSimOut);
        bandwidthVec    = cellfun(@(CELL) CELL.paramSet.bandwidth,  fetchedSimOut);
        crlbThetaVec    = cellfun(@(CELL) CELL.Jinv(1,1),           fetchedSimOut);
        %% generate MESH  
        figure;
        trisurf(delaunay(rErrVec,bandwidthVec),rErrVec,bandwidthVec,db(abs(crlbThetaVec)));   
    end
end
end