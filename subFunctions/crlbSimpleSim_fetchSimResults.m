function [fetchedSimOut] = crlbSimpleSim_fetchSimResults(cfgIn)
%% path handling
[funcpath, ~, ~]    = fileparts(mfilename('fullpath'));

%% configure
defautCfg.simOutDir             = fullfile(funcpath, '..', 'simOut', 'crlbSimpleSim');
defautCfg.filename              = 'crlbSimpleSimOut_nSensors_3_c_300000kmSec_fSample_20000Mhz.mat';
defautCfg.constantParams.names  = {'backoffFactor', 'theta' };
defautCfg.constantParams.values = {1,               pi/3    };

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

%% fetch
try
    simOut  = finalCfg.simOut;
catch
    simOutFullPath  = fullfile(finalCfg.simOutDir, finalCfg.filename);
    simOut = load(simOutFullPath);
end

try
    simOut.crlbSimpleSim_output.simOutData_CELL{1};
catch
    simOut.crlbSimpleSim_output.simOutData_CELL = cellfun(...
        @(elID) simOut.crlbSimpleSim_output.simOutData_CELL(elID), ...
        num2cell(1:numel(simOut.crlbSimpleSim_output.simOutData_CELL)), ...
        'UniformOutput', false);
end

fetchCellsIndicatorMat  = zeros(numel(finalCfg.constantParams.names),numel(simOut.crlbSimpleSim_output.simOutData_CELL));
for constPrmID = 1:numel(finalCfg.constantParams.names)
    constPrmName    = finalCfg.constantParams.names{constPrmID};
    constPrmVal     = finalCfg.constantParams.values{constPrmID};
    
    fetchCellsIndicatorMat(constPrmID,:)    = cellfun(@(CELL) CELL.paramSet.(constPrmName) == constPrmVal, simOut.crlbSimpleSim_output.simOutData_CELL);
end
fetchCellsIndicatorVec = sum(double(fetchCellsIndicatorMat)) > 0;

fetchedSimOut   = simOut.crlbSimpleSim_output.simOutData_CELL(fetchCellsIndicatorVec);

end