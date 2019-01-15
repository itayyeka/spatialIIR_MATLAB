clear all;
close all;
clc;


if true
    %% two freq approach
    if true
        %% symbolics
        syms N f DF u t r positive;
        
        f2 = f+DF;
        DU = u*(1+DF/f);
        u2 = u+DU;
        
        %% generating basic terms
        a1d1    = (r/N)*exp(1i*(u*(N-1)/2))*sin(N*u/2)/sin(u/2);
        a2d2    = (r/N)*exp(1i*(u2*(N-1)/2))*sin(N*u2/2)/sin(u2/2);
        b1d1    = (-r/N)*exp(1i*(u*(N-1)/2))*sin(N*u/2)/sin(u/2);
        b2d2    = (r/N)*exp(1i*((u2*(N-1)/2)-2*pi*DF*t))*sin(N*u2/2)/sin(u2/2);
        
        %% full expression transfer function
        hNum    = (a1d1*conj(a2d2));
        hDen    = ((1-b1d1*exp(-1i*2*pi*f*t))*(1-conj(b2d2)*exp(1i*2*pi*f2*t)));
        
        hOrig       = hNum/hDen;
        curExpr     = hOrig;
        foundExpr   = 0;
        while ~foundExpr
            [num,den]   = numden(curExpr);
            try
                num         = subs(num, u, 0);
                den         = subs(den, u, 0);
                if ~(den==0)
                    curExpr     = num/den;
                    foundExpr   = 1;
                else
                    [num,den]   = numden(curExpr);
                    num         = diff(num, u);
                    den         = diff(den, u);
                    curExpr     = num/den;
                end
            catch
                num         = diff(num, u);
                den         = diff(den, u);
                curExpr     = num/den;
            end
        end
        
        hIdeal  = subs(curExpr,t,0);
        hRel    = hOrig/hIdeal;
        hAbs    = simplify(hIdeal*conj(hIdeal));
        hRelAbs = hRel*conj(hRel);
        
        %% sanity
        if false
            curExpr     = hRelAbs;
            foundExpr   = 0;
            while ~foundExpr
                [num,den]   = numden(curExpr);
                try
                    num         = subs(num, u, 0);
                    den         = subs(den, u, 0);
                    if ~(den==0)
                        curExpr     = num/den;
                        foundExpr   = 1;
                    else
                        [num,den]   = numden(curExpr);
                        num         = diff(num, u);
                        den         = diff(den, u);
                        curExpr     = num/den;
                    end
                catch
                    num         = diff(num, u);
                    den         = diff(den, u);
                    curExpr     = num/den;
                end
            end
            
            hSanity     = subs(curExpr,t,0);
            assert(eval(hSanity)==1,'Should be 1 when in ideal scenario');
        end
        %% d/du
        curExpr         = diff(hRelAbs,u);
        foundExpr       = 0;
        while ~foundExpr
            [num,den]   = numden(curExpr);
            try
                num         = subs(num, u, 0);
                den         = subs(den, u, 0);
                if ~(den==0)
                    curExpr     = num/den;
                    foundExpr   = 1;
                else
                    [num,den]   = numden(curExpr);
                    num         = diff(num, u);
                    den         = diff(den, u);
                    curExpr     = num/den;
                end
            catch
                num         = diff(num, u);
                den         = diff(den, u);
                curExpr     = num/den;
            end
        end
        hRelAbs_d_du = subs(curExpr,t,0);
    end
    
    %% classic ULA beamwidth
    %     close all;
    %     load census;
    %     f=fit(cdate,pop,'poly2')
    %     plot(f,cdate,pop)
    
    if false
        syms p;
        syms N;
        
        orderVec= 1:1:8;
        
        for order = orderVec
            origFun         = ((sin(N*p/2))/(N*sin(p/2)))^order;
            origFun         = ((sin(N*p))/(N*sin(p)))^order;
            
            %% d_dp
            origFun_d_dp    = diff(origFun,p);
            [num,den]       = numden(origFun_d_dp);
            
            curNum      = num;
            curDen      = den;
            numDiffVal  = subs(curNum,p,0);
            denDiffVal  = subs(curDen,p,0);
            
            nIter       = 0;
            while numDiffVal == 0 && denDiffVal == 0
                curNum = diff(curNum,p);
                curDen = diff(curDen,p);
                
                numDiffVal = subs(curNum,p,0);
                denDiffVal = subs(curDen,p,0);
                
                nIter = nIter + 1;
            end
            
            lim_d_dp    = simplify((numDiffVal/denDiffVal));
            
            %% d_dp2
            origFun_d_dp2   = diff(origFun_d_dp,p);
            [num,den]       = numden(origFun_d_dp2);
            
            curNum      = num;
            curDen      = den;
            numDiffVal  = subs(curNum,p,0);
            denDiffVal  = subs(curDen,p,0);
            
            nIter       = 0;
            while numDiffVal == 0 && denDiffVal == 0
                curNum = diff(curNum,p);
                curDen = diff(curDen,p);
                
                numDiffVal = subs(curNum,p,0);
                denDiffVal = subs(curDen,p,0);
                
                nIter = nIter + 1;
            end
            
            lim_d_dp2    = simplify((numDiffVal/denDiffVal));
            
            k = simplify(lim_d_dp2/(N^2-1));
            
            kk = ((2^(-order/2)-1)/(k/2))^0.5;
            
            kkk(order) = eval(kk/pi);
            
            %% d_dp3
            origFun_d_dp3   = diff(origFun_d_dp2,p);
            [num,den]       = numden(origFun_d_dp3);
            
            curNum      = num;
            curDen      = den;
            numDiffVal  = subs(curNum,p,0);
            denDiffVal  = subs(curDen,p,0);
            
            nIter       = 0;
            while numDiffVal == 0 && denDiffVal == 0
                curNum = diff(curNum,p);
                curDen = diff(curDen,p);
                
                numDiffVal = subs(curNum,p,0);
                denDiffVal = subs(curDen,p,0);
                
                nIter = nIter + 1;
            end
            
            lim_d_dp3    = simplify((numDiffVal/denDiffVal));
            
            %% d_dp3
            origFun_d_dp4   = diff(origFun_d_dp3,p);
            [num,den]       = numden(origFun_d_dp4);
            
            curNum      = num;
            curDen      = den;
            numDiffVal  = subs(curNum,p,0);
            denDiffVal  = subs(curDen,p,0);
            
            nIter       = 0;
            while numDiffVal == 0 && denDiffVal == 0
                curNum = diff(curNum,p);
                curDen = diff(curDen,p);
                
                numDiffVal = subs(curNum,p,0);
                denDiffVal = subs(curDen,p,0);
                
                nIter = nIter + 1;
            end
            
            lim_d_dp4    = simplify((numDiffVal/denDiffVal));
            pretty(lim_d_dp4)
            
            a = (1/factorial(4))*((N^4)/5 -(2*(N^2))/3 +7/15);
            b = (1/factorial(2))*((1/3)-(N^2)/3);
            syms c;
            
            c_HPBW          = 1-(1/sqrt(2));
            nVec_HPBW       = 2:1000;
            sqrtArg         = -(4*a*c)+b^2;
            sqrtArg_HPBW    = matlabFunction(subs(sqrtArg,c,c_HPBW));
            figure; plot(sqrt(sqrtArg_HPBW(nVec_HPBW)));
            hold on;
            b_HPBW          = matlabFunction(subs(b,c,c_HPBW));
            plot(-b_HPBW(nVec_HPBW));
            plot(-b_HPBW(nVec_HPBW)-sqrt(sqrtArg_HPBW(nVec_HPBW)));
            a_HPBW          = matlabFunction(subs(2*a,c,c_HPBW));
            plot(a_HPBW(nVec_HPBW));
            y = (-b-sqrt(sqrtArg))/(2*a);
            pretty(y)
            
            
            y_HPBW = matlabFunction(subs(y,c,c_HPBW));
            
            figure;plot(nVec_HPBW.*sqrt(y_HPBW(nVec_HPBW)));
        end
        figure; plot(kkk);
        
        f = fit(orderVec(:),kkk(:),'poly2');
        
        f0 = f(0);
        
        xVec = 0:0.01:8;
        figure; plot(xVec,f(xVec));
        hold on;
        plot(orderVec,kkk,'x');
        legend({'fitted data' 'calculated factor'});
        title({...
            '3_{dB} factor' ...
            'calculated acoording to taylor expansion of various BP orders' ...
            ['reaches to ' num2str(f0) ' at 0'] ...
            'in Van-Trees - should be around 0.891' ...
            });
    end
    
    %% improvement factor
    if false
        f_improveFactor     = @(r) ((0.891)./(1-r)).*sqrt((1+2*r)/12);
        rVec                = linspace(0,1,1000);
        improveFactorVec    = f_improveFactor(rVec);
        f_rMin              = @(r) (f_improveFactor(r)-1).^2;
        rMin                = fminbnd(f_rMin, 0, 1);
        figure;
        plot(rVec, [zeros(length(improveFactorVec(:))) db(improveFactorVec(:))]);
        hold on;
        yLimVex             = ylim; % current y-axis limits
        plot([rMin rMin],[yLimVex(1) yLimVex(2)]);
        title({'Improvement factor for perfectly phase aligned feedback-based system' ...
            ['r = ' num2str(rMin) ' achieves passive ULA performance'] ...
            })
    end
    
    %% numeric validation
    if false
        f_fullExpr  = @(N,r,x,t) ...
            (1/(1-r)^2) ...
            *...
            (...
            ((N./f_D(N,x/2)).^2) - (2*r*N*cos(t+((N-1)*x/2))./f_D(N,x/2)) + r^2 ...
            );
        
        f_approxExpr    = @(N,r,x,t) ...
            1 ...
            +...
            (1/((1-r)^2))* ...
            (...
            ((1/12)*(N-1)*((1+2*r)*N-4*r+1))*(x.^2) ...
            + ...
            r*(t^2) ...
            + ...
            r*(N-1).*x.*t ...
            );
        
        nVec    = 4:8;% : 15;
        phVec   = linspace(0,.2,100);
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
    if false
        %% symbolic def
        syms x;
        syms N;
        syms r;
        syms t;
        
        origFun     = ...
            ((N*sin(x/2)/sin(N*x/2))^2) ...
            -2*N*r*cos(t+((N-1)*x/2))*sin(x/2)/sin(N*x/2) ...
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
        
        lim_d_dt2_fx = simplify((numDiffVal/denDiffVal));
        
        [num,den]       = numden(lim_d_dt2_fx);
        
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
        
        lim_d_dt2 = simplify((numDiffVal/denDiffVal));
        %% d/dxdt
        origFun_d_dxdt  = diff(origFun_d_dt,x);
        
        syms R w;
        
        origFun_d_dxdt_rw   = subs(origFun_d_dxdt,{x t}, {R*cos(w) R*sin(w)});
        
        [num,den]       = numden(origFun_d_dxdt_rw);
        
        curNum      = num;
        curDen      = den;
        numDiffVal  = subs(curNum,R,0);
        denDiffVal  = subs(curDen,R,0);
        
        nIter       = 0;
        while numDiffVal == 0 && denDiffVal == 0
            curNum      = diff(curNum,R);
            curDen      = diff(curDen,R);
            numDiffVal  = subs(curNum,R,0);
            denDiffVal  = subs(curDen,R,0);
            
            nIter = nIter + 1;
        end
        
        lim_d_dxdt = simplify((numDiffVal/denDiffVal));
    end
end

function limVal = f_D(N,x)
limVal          = sin(N*x)./sin(x);
limVal(x==0)    = N;
end