function [simOut] = beamwidthCalc(simVars, cfgIn)
simOut  = [];
%% configure
defautCfg.dummy     = 1;
defautCfg.RangeVal  = 1000;
defautCfg.cVal      = 3e8;
defautCfg.fSample   = 10e9;
lambda              = defautCfg.cVal/defautCfg.fSample;
defautCfg.DVal      = lambda/2;

cfgFields = fieldnames(defautCfg);

try
    finalCfg    = cfgIn;
catch
end

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

simOutput       = [];
simOutput.cfg   = finalCfg;
structToVariables(simVars);

%% calculate
nTheta              = 1000; 
nSensors            = length(d); 
symVarVec           = [D,       c,      omega       ];
symVarValues        = [DVal,    cVal,   2*pi*fSample];
d0                  = subs(d(:),symVarVec,symVarValues);
thetaVec            = linspace(0,2*pi,nTheta);
d0_MAT              = eval(subs(d0,theta,thetaVec));
wantedResponseMAT   = ones(nTheta)-diag(ones(nTheta,1));


noFeedbackSim       = transpose(conj(d0_MAT))*d0_MAT/nSensors;
noFeedbackSimAbs    = abs(noFeedbackSim);
% figure; plot(noFeedbackSimAbs(:,1));
% figure; surf(db(noFeedbackSimAbs));
d3dBIdxVec_noFeedback   = abs(cellfun(@(rowID) f_findFirst3Db(circshift(transpose(noFeedbackSimAbs(rowID,:)),rowID),0), num2cell(1:size(noFeedbackSimAbs,1))));
beamwidth_noFeedback    = thetaVec(d3dBIdxVec_noFeedback);

% figure; plot(beamwidth_noFeedback);
syms x;
rPole = 0.9;
feedbackSim_abs         = cell2mat(cellfun(@(doa) f_calc_iir_repsonse(fSample*DVal*cos(doa)/cVal,x,d0_MAT,rPole), num2cell(thetaVec), 'UniformOutput', false));
d3dBIdxVec_withFeedback = cellfun(@(colID) f_findFirst3Db(circshift(feedbackSim_abs(:,colID),colID-1),1), num2cell(1:size(feedbackSim_abs,1)));
beamwidth_withFeedback  = thetaVec(d3dBIdxVec_withFeedback);

figure; plot([beamwidth_noFeedback(:) beamwidth_withFeedback(:)]);
legend({'no feedback' 'with feedback'});
end

function [first3DbID] = f_findFirst3Db(responseVec,plotEnable)
responseVec_norm    = responseVec/responseVec(1);
responseVec_norm    = responseVec;
first3DbID          = find(responseVec_norm<0.5,1);
if isempty(first3DbID)
    first3DbID = 1;
end
try
    plotEnable;
catch
    plotEnable  = 0; 
end
if plotEnable %% DEBUG
    %% DEUBG
    close all;
    thetaVec    = linspace(0,2,length(responseVec_norm));
    figure; plot(thetaVec,responseVec_norm);
    disp(['last index was ' num2str(first3DbID)]);
end
end

function [response_norm] = f_calc_iir_repsonse(spatialFreq,x,d0_MAT,rPole)
order       = size(d0_MAT,1)-1;
num         = (x-exp(-1i*2*pi*spatialFreq))^order;
numCoeffs   = eval(coeffs(num,x));
den         = (x-rPole*exp(-1i*2*pi*spatialFreq))^order;
denCoeffs   = eval(coeffs(den,x));

maxCoefAbs  = max(abs([numCoeffs(:) ; denCoeffs(:)]));

numCoeffs_norm = numCoeffs/maxCoefAbs;
denCoeffs_norm = denCoeffs/maxCoefAbs;

response =  ...
    abs(...
    reshape(...
    (reshape(denCoeffs_norm,1,[])*d0_MAT) ...
    ./ ...
    (reshape(numCoeffs_norm,1,[])*d0_MAT) ...
    ,[],1) ...
    );

response_norm = response/max(abs(response));
% close all;
% figure;plot(db(abs(response_norm)));
% fvtool(numCoeffs,denCoeffs)
end