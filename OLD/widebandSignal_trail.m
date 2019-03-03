clear all;
close all;
clc;

syms r N u tau t w w1 w2 positive;

expr = ...
    (w.*exp(1i*((-pi/2)+w*(t-tau))))...
    ./...
    ((((N/r)*(sin(u/2))/(sin(N*u/2)))-exp(-1i*w*tau)).^2);

range       = 1e3;
lightSpeed  = 3e8;
fStart      = 100e6;
fStop       = fStart + 3e3;
nFreqs      = 10000; 
uVal        = 1e-6;
NVal        = 4; 
rVal        = 0.9;

tauVal      = range/lightSpeed;
exprVal     = simplify(subs(expr,{r, N, u, tau, t}, {rVal, NVal, uVal, tauVal, 0}));

freqVec     = linspace(fStart,fStop,nFreqs);
exprValVec  = eval(subs(exprVal,w,2*pi*freqVec));

exprIntegrationVec  = db(abs(cumsum(exprValVec)));
figure; plot(freqVec(:),exprIntegrationVec(:))
AAA=1;