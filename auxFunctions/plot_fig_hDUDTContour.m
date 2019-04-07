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
contourLevelList    = [-3 -40 -80 -100 -300];
if true
    %% r = 0
    hAbs_r0     = eval(subs(hAbs_symMAT,{r},{0}));
    figure;
    contourf(du_mat/pi, dPhi_mat/pi, db(abs(hAbs_r0).^2), 'LevelList', contourLevelList);
    colormap gray;
    %% r = 0.4
    hAbs_r04    = eval(subs(hAbs_symMAT,{r},{0.4}));
    hAbs_r06    = eval(subs(hAbs_symMAT,{r},{0.6}));
    hAbs_r07    = eval(subs(hAbs_symMAT,{r},{0.7}));
    hAbs_r08    = eval(subs(hAbs_symMAT,{r},{0.8}));
    figure;
    subplot(2,2,1);
    contourf(du_mat/pi, dPhi_mat/pi, db(abs(hAbs_r04).^2), 'LevelList', contourLevelList);
    %% r = 0.6
    subplot(2,2,2);
    contourf(du_mat/pi, dPhi_mat/pi, db(abs(hAbs_r06).^2), 'LevelList', contourLevelList);
    %% r = 0.7
    subplot(2,2,3);
    contourf(du_mat/pi, dPhi_mat/pi, db(abs(hAbs_r07).^2), 'LevelList', contourLevelList);
    %% r = 0.8
    subplot(2,2,4);
    contourf(du_mat/pi, dPhi_mat/pi, db(abs(hAbs_r08).^2), 'LevelList', contourLevelList);
    colormap(gray);
    %% multipeak
    rangeErrVec_multipeak           = linspace(-1.5,1.5,nPoints_sim);
    dPhiVec_multipeak               = rangeErrVec_multipeak*pi/0.5;
    [du_mat, dPhi_mat_multipeak]    = meshgrid(duVec, dPhiVec_multipeak);
    hAbs_symMAT_multipeak           = subs(hAbs2Rel,{dPhi DU N},{dPhi_mat_multipeak du_mat defaultCfg.nSensors});
    hAbs_r04_multipeak              = eval(subs(hAbs_symMAT_multipeak,{r},{0.4}));
    figure;
    contourf(du_mat/pi, dPhi_mat_multipeak/pi, db(abs(hAbs_r04_multipeak).^2), 'LevelList', contourLevelList);
    colormap gray;
end
end