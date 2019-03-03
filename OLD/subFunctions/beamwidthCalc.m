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

simOut       = [];
simOut.cfg   = finalCfg;
structToVariables(simVars);

%% calculate
nTheta              = 1000;
nSensors            = length(d);
order               = nSensors - 1;
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

simOut.beamwidth_noFeedback     = thetaVec(d3dBIdxVec_noFeedback);

% figure; plot(beamwidth_noFeedback);
syms x;
warning('off');

rPoleVec    = linspace(0.7,0.95,5);
for rPoleID = 1:numel(rPoleVec)
    feedbackSim_abs         = cellfun(@(doa) ...
        f_calc_iir_repsonse( ...
        order,              ... order,
        fSample*DVal/cVal,  ... spatialFreq,
        rPoleVec(rPoleID),  ... rPole,
        doa                 ... thetaPole
        ), num2cell(thetaVec), 'UniformOutput', false);
    warning('on');
    
    simOut.beamwidth_withFeedback(:,rPoleID)   = ...
        reshape(cellfun(@(thetaID) f_findFirst3Db(feedbackSim_abs{thetaID},1,thetaVec(thetaID)), num2cell(1:length(thetaVec))),[],1);
    
    % figure; plot(thetaVec/pi,simOut.beamwidth_withFeedback(:,rPoleID));
    
    r                   = rPoleVec(rPoleID);
    cosDelta            = (1-((1-r)^2)*2^(1/nSensors)+r^2)/(2*r);
    deltaTheta          = acos(cosDelta);
    cosThetaB           = -(deltaTheta*cVal/(2*pi*fSample*DVal)) + cos(thetaVec);
    
    simOut.beamwidth_withFeedback_analytic  = acos(cosThetaB); 
end

simOut.cfg.rPoleVec = rPoleVec;
% figure; plot([beamwidth_noFeedback(:) beamwidth_withFeedback(:)]);
% legend({'no feedback' 'with feedback'});
end

function [beamwidth] = f_findFirst3Db(responseVec,plotEnable,curTheta)
if isnumeric(responseVec)
    responseVec_norm    = responseVec/responseVec(1);
    responseVec_norm    = responseVec;
    first3DbID          = find(responseVec_norm<0.5,1);
    beamwidth           = first3DbID;
else
    thetaVec            = linspace(curTheta,2*pi+curTheta,1000);
    maxResponseAbs      = abs(responseVec(curTheta));
    f_responseAdjusted  = @(theta) -0.5+abs(responseVec(theta))/maxResponseAbs;
    resposnceValVec     = f_responseAdjusted(thetaVec);
    %     figure;plot(db(resposnceValVec));
    thetaLim            = find(resposnceValVec<0,1);
    %     figure; plot([thetaVec(thetaLim-1) thetaVec(thetaLim)],f_responseAdjusted([thetaVec(thetaLim-1) thetaVec(thetaLim)]));
    first3DbID          = fminbnd(@(theta) f_responseAdjusted(theta)^2, thetaVec(thetaLim-1), thetaVec(thetaLim));
    beamwidth           = first3DbID - curTheta;
end
if isempty(first3DbID)
    first3DbID = 1;
end
try
    plotEnable;
catch
    plotEnable  = 0;
end
if plotEnable %% DEBUG
    if isnumeric(responseVec)
        %% DEUBG
        close all;
        thetaVec    = linspace(0,2,length(responseVec_norm));
        figure; plot(thetaVec,responseVec_norm);
        disp(['last index was ' num2str(first3DbID)]);
    end
end
end

function [f_response] = f_calc_iir_repsonse(order,spatialFreq,rPole,thetaPole)
f_response = @(theta) ...
    1 ...
    ./ ...
    (...
    exp(1i*2*pi*spatialFreq*cos(theta)) ...
    - ...
    rPole*exp(1i*2*pi*spatialFreq*cos(thetaPole)) ...
    ).^order ...
    ;
end