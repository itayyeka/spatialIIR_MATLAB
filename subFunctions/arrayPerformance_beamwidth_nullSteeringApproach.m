clear all;
close all;
clc;


if true
    %% numeric validation
    if false
        f_fullExpr  = @(N,r,x,t) ...
            (1/(1-r)^2) ...
            *...
            (...
            ((N./f_D(N,x)).^2) - (2*r*N*cos(t-((N-1)*x/2))./f_D(N,x)) + r^2 ...
            );
        
        f_approxExpr    = @(N,r,x,t) 1+((1/(1-r)^2)).*((1/12)*(N-1)*((4-r)*N-3)*(x.^2) + r*(t^2));
        
        nVec    = 4:8;% : 15;
        phVec   = linspace(0,0.1,100);
        tVec    = phVec;
        rVec    = 0.9 : 0.01 : 0.99;
        
        fullExprMAT     = zeros([length(phVec) length(tVec) length(nVec) length(rVec)]);
        approxExprMAT   = zeros([length(phVec) length(tVec) length(nVec) length(rVec)]);
        
        for t   = tVec
            tID     = find(tVec    == t,1);
            for n   = nVec
                nID     = find(nVec    == n,1);
                for r   = rVec
                    rID     = find(rVec    == r,1);
                    fullExprMAT(...
                        :,   ...
                        tID,   ...
                        nID,   ...
                        rID    ...
                        ) = f_fullExpr(...
                        n,...N,
                        r,...r,...
                        phVec,...x,
                        t ...t
                        );
                    
                    approxExprMAT(...
                        :,   ...
                        tID,   ...
                        nID,   ...
                        rID    ...
                        ) = f_approxExpr(...
                        n,...N,
                        r,...r,...
                        phVec,...x,
                        t ...t
                        );
                    if false
                        close all;
                        figure;
                        plot(...
                            phVec, ...
                            fullExprMAT(...
                            :,   ...
                            tID,   ...
                            nID,   ...
                            rID    ...
                            )...
                            );
                        hold on;
                        plot(...
                            phVec, ...
                            approxExprMAT(...
                            :,   ...
                            tID,   ...
                            nID,   ...
                            rID    ...
                            )...
                            );
                    end
                end
            end
        end
        
        [phaseMAT, tauMat]  = meshgrid(phVec, tVec);
        
        for n   = nVec
            for r   = rVec
                close all;
                figure;
                surf(phaseMAT, tauMat, fullExprMAT(:,:,nVec == n,rVec == r),    'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none');
                hold on;
                surf(phaseMAT, tauMat, approxExprMAT(:,:,nVec == n,rVec == r),  'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none');
                surf(phaseMAT, tauMat, 2*ones(size(phaseMAT)),                  'FaceColor','b', 'FaceAlpha',0.5, 'EdgeColor','none');
                title(['N = ' num2str(n) ', r = ' num2str(r)]);
                xlabel('phase [RAD]');
                ylabel('tau [RAD]');
                legend({'full expression' 'approx' '3dB'});
            end
        end
    end
    %% symbolic calc
    if true
        %% symbolic def
        syms x;
        syms N;
        syms r;
        syms t;
        
        origFun     = ...
            ((N*sin(x)/sin(N*x))^2) ...
            -2*r*cos(t-((N-1)*x/2))*N*sin(x)/sin(N*x) ...
            ;
        %% d/dx
        origFun_d_dx    = diff(origFun,x);
        [num,den]       = numden(origFun_d_dx);
        
        curNum      = num;
        curDen      = den;
        numDiffVal  = subs(curNum,x,0);
        denDiffVal  = subs(curDen,x,0);
        
        nIter       = 0;
        while numDiffVal == 0 && denDiffVal == 0
            curNum = diff(curNum,x);
            curDen = diff(curDen,x);
            
            numDiffVal = subs(curNum,x,0);
            denDiffVal = subs(curDen,x,0);
            
            nIter = nIter + 1;
        end
        
        lim_d_dx    = simplify((numDiffVal/denDiffVal));
        
        %% d/dx^2
        origFun_d_dx2   = diff(origFun_d_dx,x);
        [num,den]       = numden(origFun_d_dx2);
        
        curNum      = num;
        curDen      = den;
        numDiffVal  = subs(curNum,x,0);
        denDiffVal  = subs(curDen,x,0);
        
        nIter       = 0;
        while numDiffVal == 0 && denDiffVal == 0
            curNum = diff(curNum,x);
            curDen = diff(curDen,x);
            
            numDiffVal = subs(curNum,x,0);
            denDiffVal = subs(curDen,x,0);
            
            nIter = nIter + 1;
        end
        
        lim_d_dx2   = simplify((numDiffVal/denDiffVal));
        
        %% d/dt
        origFun_d_dt    = diff(origFun,t);
        [num,den]       = numden(origFun_d_dt);
        
        curNum      = num;
        curDen      = den;
        numDiffVal  = subs(curNum,t,0);
        denDiffVal  = subs(curDen,t,0);
        
        nIter       = 0;
        while numDiffVal == 0 && denDiffVal == 0
            curNum = diff(curNum,t);
            curDen = diff(curDen,t);
            
            numDiffVal = subs(curNum,t,0);
            denDiffVal = subs(curDen,t,0);
            
            nIter = nIter + 1;
        end
        
        lim_d_dt = simplify((numDiffVal/denDiffVal));
        
        %% d/dt^2
        origFun_d_dt2   = diff(origFun_d_dt,t);
        [num,den]       = numden(origFun_d_dt2);
        
        curNum      = num;
        curDen      = den;
        numDiffVal  = subs(curNum,t,0);
        denDiffVal  = subs(curDen,t,0);
        
        nIter       = 0;
        while numDiffVal == 0 && denDiffVal == 0
            curNum = diff(curNum,t);
            curDen = diff(curDen,t);
            
            numDiffVal = subs(curNum,t,0);
            denDiffVal = subs(curDen,t,0);
            
            nIter = nIter + 1;
        end
        
        lim_d_dt2 = simplify((numDiffVal/denDiffVal));
        
        %% d/dxdt
        origFun_d_dxdt  = diff(origFun_d_dt,x);
        
        syms r w;
        
        origFun_d_dxdt_rw   = subs(origFun_d_dxdt,{x t}, {r*cos(w) r*sin(w)});
        
        [num,den]       = numden(origFun_d_dxdt_rw);
        
        curNum      = num;
        curDen      = den;
        numDiffVal  = subs(curNum,r,0);
        denDiffVal  = subs(curDen,r,0);
        
        nIter       = 0;
        while numDiffVal == 0 && denDiffVal == 0
            curNum      = diff(curNum,r);
            curDen      = diff(curDen,r);
            numDiffVal  = subs(curNum,r,0);
            denDiffVal  = subs(curDen,r,0);
            
            nIter = nIter + 1;
        end
        
        lim_d_dt2 = simplify((numDiffVal/denDiffVal));
    end
end

function limVal = f_D(N,x)
limVal          = sin(N*x)./sin(x);
limVal(x==0)    = N;
end