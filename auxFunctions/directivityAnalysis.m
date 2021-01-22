function [] = directivityAnalysis()
%% configure
nVec    = 2:10;
rVec    = [1/10 1/5 1/3 1/2 5/8 3/4 4/5];

%% symbolics
syms x real;
syms r positive;

D = @(x,N) sin(N*x)/(N*sin(x));

hNum = @(t,p,N,r) (1-r)*D(t/2,N);
hDen = @(t,p,N,r) exp(1i*(p+(N-1)*t/2)) - r*D(t/2,N);
h = @(t,p,N,r) hNum(t,p,N,r)./hDen(t,p,N,r);

%% simulate
integralResultMAT   = zeros(length(rVec),length(nVec));
exprMAT             = zeros(size(integralResultMAT));
[r_GRID, n_GRID]    = meshgrid(rVec, nVec);
for rVal = rVec
    for nVal = nVec    
        curIntRes = ...
            vpaintegral(h(x,0,nVal,rVal), 0, 2*pi...
            ) ...
            /(2*pi);
        Dir = 1 / curIntRes;
        curIntRes_eval = eval(curIntRes);
        curIntRes_STR = num2str(real(curIntRes_eval));
        Dir_eval = eval(Dir);
        Dir_STR = rats(real(Dir_eval));
        tol = abs(eval(Dir_STR)-Dir_eval);
        disp([...
            'r = ', num2str(rVal) ...
            ,', N = ', num2str(nVal) ...
            ,', int = ',  curIntRes_STR ...
            ,', Dir = ', Dir_STR ...
            ]);
        integralResultMAT(...
            rVec == rVal,...
            nVec == nVal...
            ) = real(Dir_eval);
        exprMAT(...
            rVec == rVal,...
            nVec == nVal...
            ) = (nVal-rVal)/(1-rVal);
    end    
end
figure;
surf(r_GRID, n_GRID, transpose(integralResultMAT), 'EdgeColor', "interp", 'FaceAlpha',0.75);
hold on;
plot3(r_GRID, n_GRID, transpose(exprMAT), 'dk', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'k');
legend({'Numeric integration', 'Analytic expression'});
