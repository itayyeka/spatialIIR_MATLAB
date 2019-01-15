clear all;
close all;
clc;


if true
    %% two freq approach
    if false
        %% symbolics
        syms N f DF DU t d c tS DT r positive;
        p       = pS + DP;
        pE      = d*(cos(p)-cos(pS))/c;
        Dp1     = w*pE;
        Dp2     = (w+dW)*pE;
        t       = tS + DT;
        
        %% generating the full expression
        a1d1    = (1/N)*(1-exp(1i*N*Dp1))/(1-exp(1i*Dp1));
        a2d2    = (1/N)*(1-exp(1i*N*Dp2))/(1-exp(1i*Dp2));
        b1d1    = (-1/(r*N))*(1-exp(1i*N*Dp1))/(1-exp(1i*Dp1));
        b2d2    = (1/(r*N))*((1-exp(1i*N*Dp2))/(1-exp(1i*Dp2)))*exp(-1i*dW*tS);
        
        transferFunction     = ...
            (a1d1*conj(a2d2)) ...
            / ...
            (1-(conj(b2d2)*exp(1i*(w+dW)*t))-(b1d1*exp(-1i*w*t))+b1d1*conj(b2d2)*exp(-1i*dW*t));
        
        tfAbs2  =  simplify(transferFunction*conj(transferFunction));
        
        %% tfAbs2_IDEAL
        origFun_tfAbs2_IDEAL    = tfAbs2;
        [num,den]               = numden(origFun_tfAbs2_IDEAL);
        
        curNum      = num;
        curDen      = den;
        numDiffVal  = subs(curNum,DP,0);
        denDiffVal  = subs(curDen,DP,0);
        
        nIter       = 0;
        while numDiffVal == 0 && denDiffVal == 0
            curNum = diff(curNum,DP);
            curDen = diff(curDen,DP);
            
            numDiffVal = subs(curNum,DP,0);
            denDiffVal = subs(curDen,DP,0);
            
            nIter = nIter + 1;
        end
        
        lim_tfAbs2      = (numDiffVal/denDiffVal);
        tfAbs2_IDEAL    = simplify(subs(lim_tfAbs2, DT, 0));
        
        %% tf ratio
        origFun_tfRatio     = tfAbs2_IDEAL/tfAbs2;
        [num,den]           = numden(origFun_tfRatio);
        
        curNum      = num;
        curDen      = den;
        numDiffVal  = subs(curNum,DP,0);
        denDiffVal  = subs(curDen,DP,0);
        
        nIter       = 0;
        while numDiffVal == 0 && denDiffVal == 0
            curNum = diff(curNum,DP);
            curDen = diff(curDen,DP);
            
            numDiffVal = subs(curNum,DP,0);
            denDiffVal = subs(curDen,DP,0);
            
            nIter = nIter + 1;
        end
        
        lim_tfRatio     = numDiffVal/denDiffVal;
        tfRatio_IDEAL   = simplify(subs(lim_tfRatio, DT, 0));
        
        %% d/dDP
        origFun_tfRatio_dDP = diff(origFun_tfRatio,DP);
        [num,den]           = numden(origFun_tfRatio_dDP);
        
        curNum      = num;
        curDen      = den;
        numDiffVal  = subs(curNum,DP,0);
        denDiffVal  = subs(curDen,DP,0);
        
        nIter       = 0;
        while numDiffVal == 0 && denDiffVal == 0
            curNum = diff(curNum,DP);
            curDen = diff(curDen,DP);
            
            numDiffVal = subs(curNum,DP,0);
            denDiffVal = subs(curDen,DP,0);
            
            nIter = nIter + 1;
        end
        
        lim_tfRatio_dDP     = simplify(subs(numDiffVal/denDiffVal,{DT pS},{0 pi/4}));
        
        %% d/dp2
        origFun_tfRatio_dDP2    = diff(origFun_tfRatio_dDP,DP);
        [num,den]               = numden(origFun_tfRatio_dDP2);
        
        curNum      = num;
        curDen      = den;
        numDiffVal  = subs(curNum,DP,0);
        denDiffVal  = subs(curDen,DP,0);
        
        nIter       = 0;
        while numDiffVal == 0 && denDiffVal == 0
            curNum = diff(curNum,DP);
            curDen = diff(curDen,DP);
            
            numDiffVal = subs(curNum,DP,0);
            denDiffVal = subs(curDen,DP,0);
            
            nIter = nIter + 1;
        end
        
        lim_tfRatio_dDP2    = (numDiffVal/denDiffVal);
        
        %% d/dDT
        origFun_tfRatio_dDT = diff(origFun_tfRatio,DT);
        [num,den]           = numden(origFun_tfRatio_dDT);
        
        curNum      = num;
        curDen      = den;
        numDiffVal  = subs(curNum, DT, 0);
        denDiffVal  = subs(curDen, DT, 0);
        
        nIter       = 0;
        while numDiffVal == 0 && denDiffVal == 0
            curNum = diff(curNum,DT);
            curDen = diff(curDen,DT);
            
            numDiffVal = subs(curNum, DT, 0);
            denDiffVal = subs(curDen, DT, 0);
            
            nIter = nIter + 1;
        end
        
        lim_tfRatio_dDT  = simplify(numDiffVal/denDiffVal);
        
        %% d/dDT
        origFun_tfRatio_dDT2    = diff(origFun_tfRatio,DT);
        [num,den]               = numden(origFun_tfRatio_dDT2);
        
        curNum      = num;
        curDen      = den;
        numDiffVal  = subs(curNum, DT, 0);
        denDiffVal  = subs(curDen, DT, 0);
        
        nIter       = 0;
        while numDiffVal == 0 && denDiffVal == 0
            curNum = diff(curNum,DT);
            curDen = diff(curDen,DT);
            
            numDiffVal = subs(curNum, DT, 0);
            denDiffVal = subs(curDen, DT, 0);
            
            nIter = nIter + 1;
        end
        
        lim_tfRatio_dDT2 = simplify(numDiffVal/denDiffVal);
        
        %% d/dDTdDP
        syms R;
        syms O;
        
        origFun_tfRatio_RO  = subs(origFun_tfRatio,{DP DT}, {R*cos(O) R*sin(O)});
        
        origFun_tfRatio_dDPdDT  = diff(origFun_tfRatio_RO,R);
        [num,den]               = numden(origFun_tfRatio_dDPdDT);
        
        curNum      = num;
        curDen      = den;
        numDiffVal  = subs(curNum, R, 0);
        denDiffVal  = subs(curDen, R, 0);
        
        nIter       = 0;
        while numDiffVal == 0 && denDiffVal == 0
            curNum = diff(curNum,R);
            curDen = diff(curDen,R);
            
            numDiffVal = subs(curNum, R, 0);
            denDiffVal = subs(curDen, R, 0);
            
            nIter = nIter + 1;
        end
        
        lim_tfRatio_dDPdDT = simplify(numDiffVal/denDiffVal);
    end
    
    %% classic ULA beamwidth
%     close all;
%     load census;
%     f=fit(cdate,pop,'poly2')
%     plot(f,cdate,pop)
    
    if true
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
    if true
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
    if true
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