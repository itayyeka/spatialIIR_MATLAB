function [] = plot_fig_hDUDTContour(hAbs2Rel)
simCfg.nIterations  = 1;
[defaultSimOut]     = spatialIIR_singleFreq(simCfg);
defaultCfg          = defaultSimOut.cfg;
nPoints_sim         = 100;

rangeErrVec = linspace(-0.5,0.5,nPoints_sim);
dPhiVec     = rangeErrVec*pi/0.5;
thetaVec    = linspace(0,pi,nPoints_sim);
duVec       = pi*(cos(thetaVec)-cos(defaultCfg.thetaS));

[du_mat, dPhi_mat]          = meshgrid(duVec, dPhiVec);
[dTheta_mat, rangeErr_mat]  = meshgrid(thetaVec, rangeErrVec);

syms                dPhi DU N r;
hAbs_symMAT         = subs(hAbs2Rel,{dPhi DU N},{dPhi_mat du_mat defaultCfg.nSensors});
contourLevelList    = [-3 -6 -20 -40 -80 -100 -300];
if true
    %% r = 0
    hAbs_r0     = eval(subs(hAbs_symMAT,{r},{0}));
    figure; 
    contourf(dTheta_mat/pi, dPhi_mat/pi, db(abs(hAbs_r0).^2), 'LevelList', contourLevelList);
    colormap gray;
    %% r = 0.4
    hAbs_r04    = eval(subs(hAbs_symMAT,{r},{0.4}));
    figure; 
    contourf(dTheta_mat/pi, dPhi_mat/pi, db(abs(hAbs_r04).^2), 'LevelList', contourLevelList);
    colormap gray;
    %% r = 0.6
    hAbs_r06    = eval(subs(hAbs_symMAT,{r},{0.6}));
    figure; 
    contourf(dTheta_mat/pi, dPhi_mat/pi, db(abs(hAbs_r06).^2), 'LevelList', contourLevelList);
    colormap gray;    
    %% r = 0.8
    hAbs_r08    = eval(subs(hAbs_symMAT,{r},{0.8}));
    figure; 
    contourf(dTheta_mat/pi, dPhi_mat/pi, db(abs(hAbs_r08).^2), 'LevelList', contourLevelList);
    colormap gray;
    %% r = 0.9
    hAbs_r09    = eval(subs(hAbs_symMAT,{r},{0.9}));
    figure; 
    contourf(dTheta_mat/pi, dPhi_mat/pi, db(abs(hAbs_r09).^2), 'LevelList', contourLevelList);
    colormap gray;
end
end