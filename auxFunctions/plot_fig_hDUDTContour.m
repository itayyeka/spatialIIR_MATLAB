function [] = plot_fig_hDUDTContour(hAbs2Rel)
simCfg.nIterations  = 1;
[defaultSimOut]     = spatialIIR_singleFreq(simCfg);
defaultCfg          = defaultSimOut.cfg;
nPoints_sim         = 100;

rangeErrVec = linspace(-0.5,0.5,nPoints_sim);
dPhiVec     = rangeErrVec*pi/0.5;
thetaVec    = linspace(0,pi,nPoints_sim);
duVec       = pi*(cos(thetaVec)-cos(defaultCfg.thetaS));

[dPhi_mat, du_mat]          = meshgrid(dPhiVec, duVec);
[rangeErr_mat, dTheta_mat]  = meshgrid(rangeErrVec, thetaVec);

syms                dPhi DU N r;
hAbs_symMAT         = subs(hAbs2Rel,{dPhi DU N},{dPhi_mat du_mat defaultCfg.nSensors});
contourLevelList    = [-3 -6 -20 -40 -80 -100 -300];
colormap jet;
if true
    %% r = 0
    hAbs_r0     = eval(subs(hAbs_symMAT,{r},{0}));
    figure; 
    contourf(rangeErr_mat, dTheta_mat, db(abs(hAbs_r0).^2), 'LevelList', contourLevelList);
    %% r = 0.6
    hAbs_r06    = eval(subs(hAbs_symMAT,{r},{0.6}));
    figure; 
    contourf(rangeErr_mat, dTheta_mat, db(abs(hAbs_r06).^2), 'LevelList', contourLevelList);
    %% r = 0.8
    hAbs_r08    = eval(subs(hAbs_symMAT,{r},{0.8}));
    figure; 
    contourf(rangeErr_mat, dTheta_mat, db(abs(hAbs_r08).^2), 'LevelList', contourLevelList);
    %% r = 0.95
    hAbs_r095   = eval(subs(hAbs_symMAT,{r},{0.95}));
    figure; 
    contourf(rangeErr_mat, dTheta_mat, db(abs(hAbs_r095).^2), 'LevelList', contourLevelList);
end
end