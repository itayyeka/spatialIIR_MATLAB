function [simOut] = tauErrorEffect(simVars, cfgIn)
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

simOut      = [];
simOut.cfg  = finalCfg;
structToVariables(simVars);

%% simulate
if true
    %% ideal repsonse
    syms x;
    order       = size(d,1)-1;
    thetaPole   = pi/8;
    spatialFreq = fSample*DVal*cos(thetaPole)/cVal;
    expVal      = exp(-1i*2*pi*spatialFreq);
    rPole       = 0.75;
    rZero       = rPole-0.1;
    num         = (x-rZero*expVal)^order;
    numCoeffs   = eval(coeffs(num,x));
    den         = (x-rPole*expVal)^order;
    denCoeffs   = eval(coeffs(den,x));
    
    coeffsNormFactor  = denCoeffs(end);
    
    roots(fliplr(denCoeffs))
    
    numCoeffs_norm = numCoeffs/coeffsNormFactor;
    denCoeffs_norm = denCoeffs/coeffsNormFactor;
    
    syms    tauErr;
    denCoeffs_norm_withErr = f_addTauErr(denCoeffs_norm,tauErr*2*pi*fSample);
    
    nTheta              = 100;
    symVarVec           = [D,       c,      omega       ];
    symVarValues        = [DVal,    cVal,   2*pi*fSample];
    d0                  = subs(d(:),symVarVec,symVarValues);
    thetaVec            = linspace(0,2*pi,nTheta);
    d0_MAT              = eval(subs(d0,theta,thetaVec));

    response =  ...
        abs(...
        reshape(...
        (reshape(numCoeffs_norm,1,[])*d0_MAT) ...
        ./ ...
        (reshape(denCoeffs_norm_withErr,1,[])*d0_MAT) ...
        ,[],1) ...
        );
    
    response_IDEAL  = db(abs(eval(subs(response,tauErr,0))));
    figure; plot(thetaVec/pi,response_IDEAL);
    
    %% simulate error
    nErr    = 100;
    lambda  = cVal/fSample;
    errVec  = linspace(0,.0005,nErr)*lambda/cVal;
    
    response_Err  = db(abs(eval(subs(response,tauErr,errVec))));
    
    figure;plot(thetaVec/pi,response_Err);
    legend();
    
    polesAfterErrorMat = transpose(cell2mat(cellfun(@(curTauErr) roots(fliplr(eval(subs(denCoeffs_norm_withErr,tauErr,curTauErr)))), num2cell(errVec), 'UniformOutput', false)));
    figure;
    plot(exp(1i*linspace(0,2*pi,1000)));
    hold on;
    plot(polesAfterErrorMat(1,:),'+');
    plot(polesAfterErrorMat);
    
    polesAfterErrorPhaseMat = angle(polesAfterErrorMat);
    polesAfterErrorDOAMat   = abs(polesAfterErrorMat).*exp(1i*acos(cVal*polesAfterErrorPhaseMat/(DVal*2*pi*fSample)));
    figure;
    plot(exp(1i*linspace(0,2*pi,1000)));
    hold on;
    plot(polesAfterErrorDOAMat(1,:),'+');
    plot(polesAfterErrorDOAMat);
end
end
function [denCoeffs_withErr] = f_addTauErr(denCoeffs,phaseErr)
denCoeffs_withErr       = denCoeffs;
denCoeffs_withErr(end)  = denCoeffs_withErr(end)-1;
denCoeffs_withErr       = denCoeffs_withErr*exp(1i*phaseErr);
denCoeffs_withErr(end)  = denCoeffs_withErr(end)+1;
end