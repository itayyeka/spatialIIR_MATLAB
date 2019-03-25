function [] = plot_fig_HPBW(hSym,HPBW_rootValVec)
syms DU dPhi r N;

rVec    = [0 0.1 : 0.2 : 0.9];
nVec    = 2:100;

n_r = length(rVec);
n_N = length(nVec);

resultMat   = zeros(n_N,n_r);

for nVal    = nVec
    for rVal    = rVec
        curBp           = subs(hSym,{N r dPhi},{nVal rVal 0});
        if false
            %% DEBUG
            bpValVec    = eval(subs(curBp,DU,linspace(1e-4,2*pi/nVal,100)));
            curFig      = figure;
            plot(bpValVec);
            close(curFig);
        end
        f_curBp = matlabFunction((curBp-0.5)^2);
        curHPBW = fminbnd(f_curBp,0,2*pi/nVal);
        resultMat(nVec==nVal,rVec==rVal) = curHPBW*nVal;
    end
end

if true
    %% fig_feedbackULA_HPBW_Nx_vs_N_variousR
    figure;
    %{
    DU = 2*pi*(D/lambda)*u = pi*u
    according to VanTrees pi*N*D*u/lambda -> 1.4 => Nu->2*1.4/pi
    We caclulated N*DU = N*u*pi -> 2*1.4*pi/pi = 2*1.4.
    In the plot, we show resultMat/2 to fit the known 1.4 result for r=0.
    %}
    plot(resultMat/2);
    %% fig_feedbackULA_beamwidth_limit_r_dependent
    %{
    we plot the limit as a function of r.
    It fits (1-r)(-0.4r+1.4)
    %}
    figure;
    hold on;
    rLimitVec = resultMat(end,:)/2;
    plot(rVec,rLimitVec,'-*');
    approxVal = (1-rVec).*(-0.4*rVec+1.4);
    plot(rVec,approxVal,'-diamond');
    legend('simulation','fitting');
end
end