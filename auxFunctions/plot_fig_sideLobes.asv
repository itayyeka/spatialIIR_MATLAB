function [] = plot_fig_sideLobes(hSym)
syms DU dPhi r N;

rVec    = [0 0.1 : 0.2 : 0.9];
nVec    = 2:100;

rVec    = [0 0.1 : 0.2 : 0.9];
nVec    = 2:100;

n_r = length(rVec);
n_N = length(nVec);

resultMat   = zeros(n_N,n_r);

bpSideLobe = simplify(subs(hSym,{dPhi DU},{0 3*pi/N}));

for nVal    = nVec
    for rVal    = rVec
        resultMat(nVec==nVal,rVec==rVal) = eval(subs(bpSideLobe,{N,r},{nVal,rVal}));
    end
end

figure;
plot(resultMat);
end