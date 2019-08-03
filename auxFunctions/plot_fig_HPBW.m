function [] = plot_fig_HPBW(hSym)
syms DU dPhi r N;

rVec    = [0 0.1 : 0.1 : 0.9 0.95 0.99];
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
    rIdVec  = 1:2:length(rVec(rVec<=0.9));
    plot(resultMat(:,rIdVec)/2);
    %% fig_feedbackULA_beamwidth_limit_r_dependent
    %{
    we plot the limit as a function of r.
    It fits (1-r)(-0.4r+1.4)
    %}    
    rLimitVec = resultMat(end,:)/2;
    if true
        %{
        We decided that we would like to plot B(r) = 1.4/(N\dTheta) and not
        the actual values of N\dTheta. We can then express the aperture
        improvemnt more easily.
        %}
        fig = figure;
        left_color = [0 0 1];
        right_color = [0 1 0];
        set(fig,'defaultAxesColorOrder',[left_color; right_color]);
        hold on;
        plot(rVec,rLimitVec,'-*','MarkerIndices',1:2:length(rVec));
        approxVal = (1-rVec).*(-0.4*rVec+1.4);
        plot(rVec,approxVal,'-diamond','MarkerIndices',2:2:length(rVec));
        legend('simulation','fitting');
        yyaxis right
        B           = 1.4./rLimitVec;
        B_log       = log10(B); 
        plot(rVec,B_log(:),'-.b');    
    end
end
end