function [] = directivityAnalysis()
clear;
close all;
clc;
%% configure
nVec    = 2 : 4 : 20;
rVec    = [0 0.1 : 0.1 : 0.9 0.91 : 0.01 : 0.99];
%% symbolics
syms x real;
syms r positive;

%% simulate
integralResultMAT   = zeros(length(rVec),length(nVec));
xVec                = linspace(1-10^-6,-1,1000);
xVecDiff            = diff(xVec);
phVec               = acos(xVec);
phDiffVec           = diff(phVec);
for nVal = nVec
    %% auxiliary expressions
    TN              = expand(f_TN(nVal,x));
    TN_1            = expand(f_TN(nVal-1,x));
    integrandNum    = (TN-1)*(r-1)^2;
    integrandDen    = nVal^2*(x-1)+r^2*(TN-1)+nVal*r.*(1-x-TN+TN_1);
    if false
        num = integrandNum;
        den = integrandDen;
        foundExpr   = 0;
        while ~foundExpr
            try
                numVal      = subs(num, x, 1);
                denVal      = subs(den, x, 1);
                if ~(denVal==0)
                    foundExpr   = 1;
                else
                    num         = diff(num, x);
                    den         = diff(den, x);
                end
            catch
                num         = diff(num, x);
                den         = diff(den, x);
            end
        end
        curExpr             = numVal/denVal;
        lim_integrand_1     = simplify(curExpr)
    end
    integrand = integrandNum / integrandDen;
    %% r effct on integrand - lowers the overall integration - increses directivity
    for rVal = rVec
        integrand_cur_r         = simplify(subs(integrand,r,rVal));
        integrandValues         = eval(subs(integrand_cur_r,x,xVec));
        curIntegralResult_ph    = sum(integrandValues(1:end-1).*phDiffVec/pi);
        curIntegralResult       = curIntegralResult_ph; 
        integralResultMAT(...
            rVec == rVal,...
            nVec == nVal...
            ) = curIntegralResult;
    end    
end
%% calc directivity
close all;
directivityMAT  = integralResultMAT.^(-1);
[nMAT,rMAT]     = meshgrid(nVec,rVec);
surf(rMAT,nMAT,directivityMAT);
hold on;
% [fitobject,gof,output] = fit([rMAT(:) nMAT(:)],reshape(directivityMAT./nMAT,[],1),'poly24');
% plot(fitobject)
surf(rMAT,nMAT,nMAT./((1-rMAT).^((2*nMAT)./(2*nMAT-1))));

end
function [TN] = f_TN(N,x)
TN = 0;
for k = 0 : floor(N/2)
    TN = TN + nchoosek(N,2*k)*((x^2-1).^k).*x.^(N-2*k);
end
end