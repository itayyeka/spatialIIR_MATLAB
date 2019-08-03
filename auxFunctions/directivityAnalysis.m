function [] = directivityAnalysis(cfgIn)
try
    cfgIn;
catch
    clear;
    close all;
    clc;
end
%% configure
nVec    = 2:6;
rVec    = [0.25 : 0.05 : 0.75];
rVec    = [0.2 : 0.1 : 0.7];
%% symbolics
syms x real;
syms r positive;

D = @(x,N) sin(N*x)/(N*sin(x));

hNum = @(t,p,N,r) (1-r)*D(t/2,N);
hDen = @(t,p,N,r) exp(1i*(p+(N-1)*t/2)) - r*D(t/2,N);
h = @(t,p,N,r) hNum(t,p,N,r)./hDen(t,p,N,r);

%% simulate
integralResultMAT   = zeros(length(rVec),length(nVec));
for nVal = nVec
    %% r effct on integrand - lowers the overall integration - increses directivity
    for rVal = rVec
        curIntRes = ...
            vpaintegral(h(x,0,nVal,rVal), 0, 2*pi...
            ...,'IgnoreAnalyticConstraints', true, 'IgnoreSpecialCases', true, 'PrincipalValue',true ...
            ) ...
            /(2*pi);
        Dir = 1 / curIntRes;
        curIntRes_eval = eval(curIntRes);
        curIntRes_STR = num2str(real(curIntRes_eval));
        Dir_eval = eval(Dir);
        Dir_STR = rats(real(Dir_eval));
        tol = abs(eval(Dir_STR)-Dir_eval);
        expr = (rVal^2-(nVal-1)*rVal+nVal)/((1-rVal)^2);
        disp([...
            'r = ', num2str(rVal) ...
            ,', N = ', num2str(nVal) ...
            ,', int = ',  curIntRes_STR ...
            ,', Dir = ', Dir_STR ...
            ...,', Dir(1-r)^2 = ', rats(eval(Dir_STR)*(1-rVal)^2) ...
            ...,', tol = ', num2str(tol) ...
            ', Expr = ', rats(expr), ...
            ]);
        integralResultMAT(...
            rVec == rVal,...
            nVec == nVal...
            ) = real(Dir_eval);
    end
end
integralResultMAT_CELL = num2cell(integralResultMAT);
integralResultMAT_CELL_rat = cellfun(@(x) rats(x), integralResultMAT_CELL, 'UniformOutput', false);
rVec_CELL = num2cell(rVec(:));
rVec_CELL_STR = cellfun(@(x) num2str(x), rVec_CELL, 'UniformOutput', false);
nVec_CELL_STR = cellfun(@(x) ['N' num2str(x)], num2cell(nVec), 'UniformOutput', false);
T = array2table([rVec_CELL_STR(:) integralResultMAT_CELL_rat],...
    'VariableNames',[{'r'} reshape(nVec_CELL_STR,1,[])])
