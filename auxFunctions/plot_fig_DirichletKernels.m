function [] = plot_fig_DirichletKernels()
    clear all;
    close all;
    clc;
    
    %% cfg
    nVec = [2, 5, 7, 10];
    nTheta = 100;
    
    %% aux
    f_D = @(N,x) sin(N*x)./sin(x);
    
    %% sim
    thetaVec = linspace(-pi,pi,nTheta);
    res = zeros(nTheta,numel(nVec));
    resNorm = zeros(nTheta,numel(nVec));
    for n = nVec
        res(:, nVec==n) = f_D(n,thetaVec);
        resNorm(:, nVec==n) = f_D(n,thetaVec)/n;
    end
    
    %% plot
    figure;
    plot(thetaVec/pi,res);
    axis tight;
    figure;
    plot(thetaVec/pi,resNorm);
    axis tight;
    fixfig;
end

