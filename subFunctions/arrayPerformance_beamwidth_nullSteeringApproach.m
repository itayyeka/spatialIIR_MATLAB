clear all;
close all;
clc;


if true
    %% two freq approach
    if false
        %% symbolics
        syms N f DF u t r positive;
        
        f2 = f+DF;
        DU = u*(1+DF/f);
        u2 = u+DU;
        
        %% generating basic terms
        if true
            a1d1    = (r/N)*exp(1i*(u*(N-1)/2))*sin(N*u/2)/sin(u/2);
            a2d2    = (r/N)*exp(1i*(u2*(N-1)/2))*sin(N*u2/2)/sin(u2/2);
            b1d1    = (-r/N)*exp(1i*(u*(N-1)/2))*sin(N*u/2)/sin(u/2);
            b2d2    = (r/N)*exp(1i*((u2*(N-1)/2)-2*pi*DF*t))*sin(N*u2/2)/sin(u2/2);
        else
            a1d1    = (r/N)*(1-exp(1i*N*u))/(1-exp(1i*u));
            a2d2    = (r/N)*(1-exp(1i*N*u2))/(1-exp(1i*u2));
            b1d1    = (-r/N)*(1-exp(1i*N*u))/(1-exp(1i*u));
            b2d2    = (r/N)*exp(-1i*2*pi*DF*t)*(1-exp(1i*N*u2))/(1-exp(1i*u2));
        end
        %% full expression transfer function
        hNum    = (a1d1*conj(a2d2));
        hDen    = ((1-b1d1*exp(-1i*2*pi*f*t))*(1-conj(b2d2)*exp(1i*2*pi*f2*t)));
        
        hOrig       = simplify(hNum/hDen);
        
        curExpr     = subs(hOrig, t, 0);
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
        
        hIdeal      = simplify(curExpr);
        hIdealAbs   = simplify(hIdeal*conj(hIdeal));
        hRel        = hOrig/hIdeal;
        hRelAbs     = hRel*conj(hRel);
        hOrigAbs    = hOrig*conj(hOrig);
        
        %% sanity
        if false
            curExpr     = simplify(subs(hRelAbs, t, 0));
            foundExpr   = 0;
            [num,den]   = numden(curExpr);
            while ~foundExpr
                try
                    numVal  = subs(num, u, 0);
                    denVal  = subs(den, u, 0);
                    if ~(denVal==0)
                        foundExpr   = 1;
                    else
                        num         = diff(num, u);
                        den         = diff(den, u);
                    end
                catch
                    num         = diff(num, u);
                    den         = diff(den, u);
                end
            end
            curExpr     = numVal/denVal;
            lim_hRelAbs = simplify(curExpr)
            assert(eval(lim_hRelAbs)==1,'Should be 1 when in ideal scenario');
        end
        %% d/du
        hOrigAbs_d_du   = simplify(subs(diff(hOrigAbs,u),t,0));
        if false
            curExpr     = hOrigAbs_d_du;
            foundExpr   = 0;
            [num,den]   = numden(curExpr);
            while ~foundExpr
                try
                    numVal  = subs(num, u, 0);
                    denVal  = subs(den, u, 0);
                    if ~(denVal==0)
                        foundExpr   = 1;
                    else
                        num         = diff(num, u);
                        den         = diff(den, u);
                    end
                catch
                    num         = diff(num, u);
                    den         = diff(den, u);
                end
            end
            curExpr             = simplify(numVal/denVal);
            lim_hOrigAbs_d_du   = curExpr
        end
        %% d/du2
        hOrigAbs_d_du2  = simplify(subs(diff(hOrigAbs_d_du,u),t,0));
        if false
            curExpr     = hOrigAbs_d_du2;
            foundExpr   = 0;
            [num,den]   = numden(curExpr);
            while ~foundExpr
                try
                    numVal  = subs(num, u, 0);
                    denVal  = subs(den, u, 0);
                    if ~(denVal==0)
                        foundExpr   = 1;
                    else
                        num         = diff(num, u);
                        den         = diff(den, u);
                    end
                catch
                    num         = diff(num, u);
                    den         = diff(den, u);
                end
            end
            curExpr             = simplify(numVal/denVal);
            lim_hOrigAbs_d_du2  = curExpr
        end
        %% d/dt
        hOrigAbs_d_dt   = simplify(subs(diff(hOrigAbs,t),t,0));
        if false
            curExpr     = hOrigAbs_d_dt;
            foundExpr   = 0;
            [num,den]   = numden(curExpr);
            while ~foundExpr
                try
                    numVal  = subs(num, u, 0);
                    denVal  = subs(den, u, 0);
                    if ~(denVal==0)
                        foundExpr   = 1;
                    else
                        num         = diff(num, u);
                        den         = diff(den, u);
                    end
                catch
                    num         = diff(num, u);
                    den         = diff(den, u);
                end
            end
            curExpr             = simplify(numVal/denVal);
            lim_hOrigAbs_d_dt   = curExpr
        end
        %% d/dt2
        hOrigAbs_d_dt2  = simplify(subs(diff(hOrigAbs_d_dt,t),t,0));
        if false
            curExpr     = hOrigAbs_d_dt2;
            foundExpr   = 0;
            [num,den]   = numden(curExpr);
            while ~foundExpr
                try
                    numVal  = subs(num, u, 0);
                    denVal  = subs(den, u, 0);
                    if ~(denVal==0)
                        foundExpr   = 1;
                    else
                        num         = diff(num, u);
                        den         = diff(den, u);
                    end
                catch
                    num         = diff(num, u);
                    den         = diff(den, u);
                end
            end
            curExpr             = simplify(numVal/denVal);
            lim_hOrigAbs_d_dt2  = curExpr
        end
        %% d/dt2
        syms R w positive;
        hOrigAbs_d_dudt  = simplify(subs(diff(hOrigAbs_d_dt,u),{u t},{R*cos(w) R*sin(w)}));
        if false
            curExpr     = hOrigAbs_d_dudt;
            foundExpr   = 0;
            [num,den]   = numden(curExpr);
            while ~foundExpr
                try
                    numVal  = subs(num, R, 0);
                    denVal  = subs(den, R, 0);
                    if ~(denVal==0)
                        foundExpr   = 1;
                    else
                        num         = diff(num, R);
                        den         = diff(den, R);
                    end
                catch
                    num         = diff(num, R);
                    den         = diff(den, R);
                end
            end
            curExpr             = simplify(numVal/denVal);
            lim_hOrigAbs_d_dudt = curExpr
        end
    end
    
    %% classic ULA beamwidth
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
    
    %% symbolic calc - no attenuation
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
    
    %% numeric validation  - no attenuation
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
    
    %% symbolic calc - with attenuation
    if true
        %% symbolics
        syms N DU ts DT c w r positive;
        syms k;
        
        t   = ts + DT;
        gs  = 1/(c*ts);
        g   = 1/(c*t);
        
        %% generating basic terms
        ad  = (r/N)*exp(1i*(DU*(N-1)/2)).*sin(N*DU/2)./sin(DU/2);
        bd  = ad * exp(1i*w*ts) / (gs^2) ;
        
        %% full expression transfer function
        h       = (ad*(gs^2))/(1-(g^2)*bd*exp(-1i*w*t));
        hAbs2   = simplify(h*conj(h));
        
        %% Ideal
        if true
            curExpr     = subs(hAbs2, DT, 0);
            [num,den]   = numden(curExpr);
            foundExpr   = 0;
            while ~foundExpr
                try
                    numVal      = subs(num, DU, 0);
                    denVal      = subs(den, DU, 0);
                    if ~(denVal==0)
                        foundExpr   = 1;
                    else
                        num         = simplify(diff(num, DU));
                        den         = simplify(diff(den, DU));
                    end
                catch
                    num         = simplify(diff(num, DU));
                    den         = simplify(diff(den, DU));
                end
            end
            curExpr     = numVal/denVal;
            lim_hAbs2   = simplify(curExpr)
        end
        hAbs2Rel        = lim_hAbs2/hAbs2;
        
        %% taylor
        if false
            syms R phi;
            %% d/dDU : 0
            hAbs2Rel_dDU    = diff(hAbs2Rel,DU);
            if false
                curExpr     = subs(hAbs2Rel_dDU,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, DU));
                            den         = simplify(diff(den, DU));
                        end
                    catch
                        num         = simplify(diff(num, DU));
                        den         = simplify(diff(den, DU));
                    end
                end
                curExpr         = numVal/denVal;
                lim_hAbs2Rel_du = simplify(curExpr)
            end
            
            %% d/dDT : 4*r/(ts*(r-1))
            hAbs2Rel_dDT    = diff(hAbs2Rel,DT);
            if false
                curExpr     = subs(hAbs2Rel_dDT,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, DU));
                            den         = simplify(diff(den, DU));
                        end
                    catch
                        num         = simplify(diff(num, DU));
                        den         = simplify(diff(den, DU));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDT    = simplify(curExpr)
            end
            
            %% d/dDU2 : ((N-1)*(N-4*r+2*N*r*+1))/(6*(r-1)^2)
            hAbs2Rel_dDU2   = diff(hAbs2Rel_dDU,DU);
            if false
                curExpr     = subs(hAbs2Rel_dDU2,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, DU));
                            den         = simplify(diff(den, DU));
                        end
                    catch
                        num         = simplify(diff(num, DU));
                        den         = simplify(diff(den, DU));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU2   = simplify(curExpr)
            end
            
            %% d/dDUDT : -(r*w*(N - 1))/(r - 1)^2
            hAbs2Rel_dDUDT  = diff(hAbs2Rel_dDT,DU);
            if false
                curExpr     = subs(hAbs2Rel_dDUDT,{DU DT},{R*cos(phi),R*sin(phi)});
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, R, 0);
                        denVal      = subs(den, R, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, R));
                            den         = simplify(diff(den, R));
                        end
                    catch
                        num         = simplify(diff(num, R));
                        den         = simplify(diff(den, R));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDUDT  = simplify(curExpr)
            end
            
            %% d/dDT2 : (r*(2*(tsVal*wVal)^2+20*r-12))/((tsVal*(r-1))^2)
            hAbs2Rel_dDT2 = diff(hAbs2Rel_dDT,DT);
            if false
                curExpr     = subs(hAbs2Rel_dDT2,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, DU));
                            den         = simplify(diff(den, DU));
                        end
                    catch
                        num         = simplify(diff(num, DU));
                        den         = simplify(diff(den, DU));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDT2   = simplify(curExpr)
            end
            
            %% d/dDU3 : 0
            hAbs2Rel_dDU3   = diff(hAbs2Rel_dDU2,DU);
            if false
                curExpr     = subs(hAbs2Rel_dDU3,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, DU));
                            den         = simplify(diff(den, DU));
                        end
                    catch
                        num         = simplify(diff(num, DU));
                        den         = simplify(diff(den, DU));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU3   = simplify(curExpr)
            end
            
            %% d/dDU2DT : -(2*r*(N^2 - 3*N + 2))/(3*ts*(r - 1)^2)
            hAbs2Rel_dDU2DT = diff(hAbs2Rel_dDU2,DT);
            if false
                curExpr     = subs(hAbs2Rel_dDU2DT,{DU DT},{R*cos(phi),R*sin(phi)});
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, R, 0);
                        denVal      = subs(den, R, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, R));
                            den         = simplify(diff(den, R));
                        end
                    catch
                        num         = simplify(diff(num, R));
                        den         = simplify(diff(den, R));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU2DT  = simplify(curExpr)
            end
            
            %% d/dDUDT2 : (4*r*w*(N - 1))/(ts*(r - 1)^2)
            hAbs2Rel_dDUDT2 = diff(hAbs2Rel_dDT2,DU);
            if false
                curExpr     = subs(hAbs2Rel_dDUDT2,{DU DT},{R*cos(phi),R*sin(phi)});
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, R, 0);
                        denVal      = subs(den, R, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, R));
                            den         = simplify(diff(den, R));
                        end
                    catch
                        num         = simplify(diff(num, R));
                        den         = simplify(diff(den, R));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDUDT2  = simplify(curExpr)
            end
            
            %% d/dDT3 : -(12*r*(ts^2*w^2 + 10*r - 4))/(ts^3*(r - 1)^2)
            hAbs2Rel_dDT3   = diff(hAbs2Rel_dDT2,DT);
            if false
                curExpr     = subs(hAbs2Rel_dDT3,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, DU));
                            den         = simplify(diff(den, DU));
                        end
                    catch
                        num         = simplify(diff(num, DU));
                        den         = simplify(diff(den, DU));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDT3   = simplify(expand(curExpr))
            end
            
            %% d/dDU4 : (30*n*r-12*r-20*n^2*r+2*n^4*r-5*n^2+3*n^4+2)/(30*(r-1)^2)
            hAbs2Rel_dDU4   = diff(hAbs2Rel_dDU3,DU);
            if false
                curExpr     = subs(hAbs2Rel_dDU4,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, DU));
                            den         = simplify(diff(den, DU));
                        end
                    catch
                        num         = simplify(diff(num, DU));
                        den         = simplify(diff(den, DU));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU4   = simplify(curExpr)
            end
            
            %% d/dDT3DU : (r*w*(ts^2*w^2 - 18)*(N - 1))/(ts^2*(r - 1)^2)
            hAbs2Rel_dDT3DU = diff(hAbs2Rel_dDT3,DU);
            if false
                curExpr     = subs(hAbs2Rel_dDT3DU,{DU DT},{R*cos(phi),R*sin(phi)});
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, R, 0);
                        denVal      = subs(den, R, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, R));
                            den         = simplify(diff(den, R));
                        end
                    catch
                        num         = simplify(diff(num, R));
                        den         = simplify(diff(den, R));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDT3DU  = simplify(curExpr)
            end
            
            %% d/dDU2DT2 : -(r*(ts^2*w^2 - 6)*(N^2 - 3*N + 2))/(3*ts^2*(r - 1)^2)
            hAbs2Rel_dDU2dDT2  = diff(diff(hAbs2Rel_dDUDT,DU),DT);
            if false
                curExpr     = subs(hAbs2Rel_dDU2dDT2,{DU DT},{R*cos(phi),R*sin(phi)});
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, R, 0);
                        denVal      = subs(den, R, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, R));
                            den         = simplify(diff(den, R));
                        end
                    catch
                        num         = simplify(diff(num, R));
                        den         = simplify(diff(den, R));
                    end
                end
                curExpr                 = numVal/denVal;
                lim_hAbs2Rel_dDU2dDT2   = simplify(curExpr)
            end
            
            %% d/dDTDU3 : -(r*w*(N - 1)^2)/(2*(r - 1)^2)
            hAbs2Rel_dDTDU3 = diff(hAbs2Rel_dDU3,DT);
            if false
                curExpr     = subs(hAbs2Rel_dDTDU3,{DU DT},{R*cos(phi),R*sin(phi)});
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, R, 0);
                        denVal      = subs(den, R, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, R));
                            den         = simplify(diff(den, R));
                        end
                    catch
                        num         = simplify(diff(num, R));
                        den         = simplify(diff(den, R));
                    end
                end
                curExpr                 = numVal/denVal;
                lim_hAbs2Rel_dDTDU3     = simplify(curExpr)
            end
            
            %% d/dDT4 : (r*(- 2*ts^4*w^4 + 72*ts^2*w^2 + 840*r - 240))/(ts^4*(r - 1)^2)
            hAbs2Rel_dDT4   = diff(hAbs2Rel_dDT3,DT);
            if false
                curExpr     = subs(hAbs2Rel_dDT4,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, DU));
                            den         = simplify(diff(den, DU));
                        end
                    catch
                        num         = simplify(diff(num, DU));
                        den         = simplify(diff(den, DU));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDT4   = simplify(curExpr)
            end
            
            %% d/dDU5 : 0
            hAbs2Rel_dDU5   = diff(hAbs2Rel_dDU4,DU);
            if false
                curExpr     = subs(hAbs2Rel_dDU5,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, DU));
                            den         = simplify(diff(den, DU));
                        end
                    catch
                        num         = simplify(diff(num, DU));
                        den         = simplify(diff(den, DU));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU5   = simplify(curExpr)
            end
            
            %% d/dDU4DT -(2*r*(N^4 - 10*N^2 + 15*N - 6))/(15*ts*(r - 1)^2): 
            hAbs2Rel_dDU4DT = diff(hAbs2Rel_dDU4,DT);
            if false
                curExpr     = subs(hAbs2Rel_dDU4DT,{DU DT},{R*cos(phi),R*sin(phi)});
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, R, 0);
                        denVal      = subs(den, R, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, R));
                            den         = simplify(diff(den, R));
                        end
                    catch
                        num         = simplify(diff(num, R));
                        den         = simplify(diff(den, R));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU4DT  = simplify(curExpr)
            end
            
            %% d/dDU3DT2 : (2*r*w*(N - 1)^2)/(ts*(r - 1)^2)
            hAbs2Rel_dDU3DT2 = diff(hAbs2Rel_dDU2dDT2,DU);
            if false
                curExpr     = subs(hAbs2Rel_dDU3DT2,{DU DT},{R*cos(phi),R*sin(phi)});
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, R, 0);
                        denVal      = subs(den, R, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, R));
                            den         = simplify(diff(den, R));
                        end
                    catch
                        num         = simplify(diff(num, R));
                        den         = simplify(diff(den, R));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU3DT2  = simplify(curExpr)
            end
            
            %% d/dDU2DT3 : (2*r*(ts^2*w^2 - 4)*(N^2 - 3*N + 2))/(ts^3*(r - 1)^2)
            hAbs2Rel_dDU2DT3 = diff(hAbs2Rel_dDU2dDT2,DT);
            if false
                curExpr     = subs(hAbs2Rel_dDU2DT3,{DU DT},{R*cos(phi),R*sin(phi)});
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, R, 0);
                        denVal      = subs(den, R, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, R));
                            den         = simplify(diff(den, R));
                        end
                    catch
                        num         = simplify(diff(num, R));
                        den         = simplify(diff(den, R));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU2DT3  = simplify(curExpr)
            end
            
            %% d/dDUDT4 : -(8*r*w*(ts^2*w^2 - 12)*(N - 1))/(ts^3*(r - 1)^2)
            hAbs2Rel_dDUDT4 = diff(hAbs2Rel_dDT3DU,DT);
            if false
                curExpr     = subs(hAbs2Rel_dDUDT4,{DU DT},{R*cos(phi),R*sin(phi)});
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, R, 0);
                        denVal      = subs(den, R, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, R));
                            den         = simplify(diff(den, R));
                        end
                    catch
                        num         = simplify(diff(num, R));
                        den         = simplify(diff(den, R));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDUDT4  = simplify(curExpr)
            end
            
            %% d/dDT5 : -(20*r*(- ts^4*w^4 + 24*ts^2*w^2 + 336*r - 72))/(ts^5*(r - 1)^2)
            hAbs2Rel_dDT5   = diff(hAbs2Rel_dDT4,DT);
            if true
                curExpr     = subs(hAbs2Rel_dDT5,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, DU));
                            den         = simplify(diff(den, DU));
                        end
                    catch
                        num         = simplify(diff(num, DU));
                        den         = simplify(diff(den, DU));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDT5   = simplify(expand(curExpr))
            end
            
            %% d/dDU6 : (24*r - 84*n*r + 84*n^2*r - 28*n^4*r + 4*n^6*r + 14*n^2 - 21*n^4 + 10*n^6 - 3)/(84*(r - 1)^2)
            hAbs2Rel_dDU6   = diff(hAbs2Rel_dDU5,DU);
            if false
                curExpr     = subs(hAbs2Rel_dDU6,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = simplify(diff(num, DU));
                            den         = simplify(diff(den, DU));
                        end
                    catch
                        num         = simplify(diff(num, DU));
                        den         = simplify(diff(den, DU));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU6   = simplify(curExpr)
            end
            
        end
        
        %% numeric validation
        f_hAbs2Rel  = matlabFunction(hAbs2Rel,'Vars',{[DT, N, c, r, ts, w]' DU});
        if true
            %% configure
            cVal        = 3e8;
            fVal        = 10e9;
            RVal        = 100e3;
            lambdaVal   = cVal/fVal;
            if true
                %% N
                init_N      = 3;
                final_N     = 10;
                hop_N       = 1;
                %% r
                nVal_r      = 30;
                init_r      = 0.6;
                final_r     = 0.99;
                %% DU
                nVal_DU     = 30;
                init_DU     = 1e-6;
                final_DU    = .5;
                %% DT
                nVal_DT     = 30;
                init_DT     = 0;
                final_DT    = .04*lambdaVal/cVal;
            end
            
            %% generate input
            nVec    = init_N : hop_N : final_N;
            rVec    = linspace(init_r,final_r,nVal_r);
            DUVec   = linspace(init_DU,final_DU,nVal_DU);
            DTVec   = linspace(init_DT,final_DT,nVal_DT);
            
            %% simulate
            nVal_N              = length(nVec);
            simResult           = zeros(nVal_DT,nVal_DU,nVal_N,nVal_r);
            approxResult        = zeros(nVal_DT,nVal_DU,nVal_N,nVal_r);
            approxResult_noTs   = zeros(nVal_DT,nVal_DU,nVal_N,nVal_r);
            rootMat             = zeros(nVal_N,nVal_r);
            symVarValues        = zeros(6,1);
            
            %[DT, N, c, r, ts, w]
            wVal    = 2*pi*fVal;
            tsVal   = RVal/cVal;
            symVarValues(6) = wVal;
            symVarValues(3) = cVal;
            symVarValues(5) = tsVal;
            
            for r = rVec
                symVarValues(4) = r;
                for n = nVec
                    symVarValues(2) = n;
                    cDU             = (nchoosek(1,0)/factorial(1))*(0);
                    cDT             = (nchoosek(1,1)/factorial(1))*(4*r/(tsVal*(r-1)));
                    cDT_noTs        = 0;
                    cDU2            = (nchoosek(2,0)/factorial(2))*((n-1)*(n-4*r+2*n*r*+1))/(6*(r-1)^2);
                    cDUDT           = (nchoosek(2,1)/factorial(2))*(-(r*wVal*(n - 1))/(r - 1)^2);
                    cDT2            = (nchoosek(2,2)/factorial(2))*((r*(2*(tsVal*wVal)^2+20*r-12))/((tsVal*(r-1))^2));
                    cDT2_noTs       = (nchoosek(2,2)/factorial(2))*((r*(2*(      wVal)^2        ))/((      (r-1))^2));
                    cDU3            = (nchoosek(3,0)/factorial(3))*(0);
                    cDU2DT          = (nchoosek(3,1)/factorial(3))*(-(2*r*(n^2 - 3*n + 2))/(3*tsVal*(r - 1)^2));
                    cDU2DT_noTs     = (nchoosek(3,1)/factorial(3))*(0);
                    cDUDT2          = (nchoosek(3,2)/factorial(3))*((4*r*wVal*(n - 1))/(tsVal*(r - 1)^2));
                    cDUDT2_noTs     = (nchoosek(3,2)/factorial(3))*(0);
                    cDT3            = (nchoosek(3,3)/factorial(3))*(-(12*r*(tsVal^2*wVal^2 + 10*r - 4))/(tsVal^3*(r - 1)^2));
                    cDT3_noTs       = (nchoosek(3,3)/factorial(3))*(0);
                    cDU4            = (nchoosek(4,0)/factorial(4))*(30*n*r-12*r-20*n^2*r+2*n^4*r-5*n^2+3*n^4+2)/(30*(r-1)^2);
                    cDU3DT          = (nchoosek(4,1)/factorial(4))*(-(r*wVal*(n - 1)^2)/(2*(r - 1)^2));
                    cDU2DT2         = (nchoosek(4,2)/factorial(4))*(-(r*(tsVal^2*wVal^2 - 6)*(n^2 - 3*n + 2))/(3*tsVal^2*(r - 1)^2));
                    cDU2DT2_noTs    = (nchoosek(4,2)/factorial(4))*(-(r*(        wVal^2    )*(n^2 - 3*n + 2))/(3*     (r - 1)^2));
                    cDUDT3          = (nchoosek(4,3)/factorial(4))*((r*wVal*(tsVal^2*wVal^2 - 18)*(n - 1))/(tsVal^2*(r - 1)^2));
                    cDUDT3_noTs     = (nchoosek(4,3)/factorial(4))*((r*wVal*(        wVal^2     )*(n - 1))/(     (r - 1)^2));
                    cDT4            = (nchoosek(4,4)/factorial(4))*((r*(- 2*tsVal^4*wVal^4 + 72*tsVal^2*wVal^2 + 840*r - 240))/(tsVal^4*(r - 1)^2));
                    cDT4_noTs       = (nchoosek(4,4)/factorial(4))*((r*(- 2*        wVal^4                                  ))/(        (r - 1)^2));
                    cDU5            = (nchoosek(0,5)/factorial(5))*(0);
                    cDU4DT          = (nchoosek(1,5)/factorial(5))*(0);
                    cDU3DT2         = (nchoosek(2,5)/factorial(5))*(0);
                    cDU2DT3         = (nchoosek(3,5)/factorial(5))*(0);
                    cDUDT4          = (nchoosek(4,5)/factorial(5))*(0);
                    cDT5            = (nchoosek(5,5)/factorial(5))*(0);
                    cDU6            = (nchoosek(6,0)/factorial(6))*(24*r - 84*n*r + 84*n^2*r - 28*n^4*r + 4*n^6*r + 14*n^2 - 21*n^4 + 10*n^6 - 3)/(84*(r - 1)^2);
                    cVec6           = [cDU6 cDU4 cDU2 -1]; % compare to 2
                    curRootVec6     = roots(cVec6);
                    [~, minImagID]  = min(abs(imag(curRootVec6)));
                    curRoot6        = sqrt(real(curRootVec6(minImagID)));
                    cVec4           = [cDU4 cDU2 -1]; % compare to 2
                    curRootVec4     = roots(cVec4);
                    curRoot4        = sqrt(curRootVec4(curRootVec4>0));
                    curRoot         = curRoot4;
                    rootMat(...
                        nVec    == n,       ... n
                        rVec    == r        ... r
                        ) = curRoot;
                    for curDT = DTVec
                        symVarValues(1) = curDT;
                        simResult(...
                            DTVec   == curDT,   ... DT
                            :,                  ... DU
                            nVec    == n,       ... n
                            rVec    == r        ... r
                            ) =                 ...
                            f_hAbs2Rel(symVarValues,DUVec);
                        approxResult(...
                            DTVec   == curDT,   ... DT
                            :,                  ... DU
                            nVec    == n,       ... n
                            rVec    == r        ... r
                            ) =                 ...
                            1 ...
                            +cDU            *(DUVec.^(1))*curDT^(0) ...
                            +cDT            *(DUVec.^(0))*curDT^(1) ...
                            ...+cDT_noTs       *(DUVec.^(0))*curDT^(1) ...
                            +cDU2           *(DUVec.^(2))*curDT^(0) ...
                            +cDUDT          *(DUVec.^(1))*curDT^(1) ...
                            +cDT2           *(DUVec.^(0))*curDT^(2) ...
                            ...+cDT2_noTs      *(DUVec.^(0))*curDT^(2) ...
                            +cDU3           *(DUVec.^(3))*curDT^(0) ...
                            +cDU2DT         *(DUVec.^(2))*curDT^(1) ...
                            ...+cDU2DT_noTs    *(DUVec.^(2))*curDT^(1) ...
                            +cDUDT2         *(DUVec.^(1))*curDT^(2) ...
                            ...+cDUDT2_noTs    *(DUVec.^(1))*curDT^(2) ...
                            +cDT3           *(DUVec.^(0))*curDT^(3) ...
                            ...+cDT3_noTs      *(DUVec.^(0))*curDT^(3) ...
                            +cDU4           *(DUVec.^(4))*curDT^(0) ...
                            +cDU3DT         *(DUVec.^(3))*curDT^(1) ...
                            +cDU2DT2        *(DUVec.^(2))*curDT^(2) ...
                            ...+cDU2DT2_noTs   *(DUVec.^(2))*curDT^(2) ...
                            +cDUDT3         *(DUVec.^(1))*curDT^(3) ...
                            ...+cDUDT3_noTs    *(DUVec.^(1))*curDT^(3) ...
                            +cDT4           *(DUVec.^(0))*curDT^(4) ...
                            ...+cDT4_noTs      *(DUVec.^(0))*curDT^(4) ...
                            ...+cDU6           *(DUVec.^(6))*curDT^(0) ...
                            ;
                        approxResult_noTs(...
                            DTVec   == curDT,   ... DT
                            :,                  ... DU
                            nVec    == n,       ... n
                            rVec    == r        ... r
                            ) =                 ...
                            1 ...
                            +cDU            *(DUVec.^(1))*curDT^(0) ...
                            ...+cDT            *(DUVec.^(0))*curDT^(1) ...
                            +cDT_noTs       *(DUVec.^(0))*curDT^(1) ...
                            +cDU2           *(DUVec.^(2))*curDT^(0) ...
                            +cDUDT          *(DUVec.^(1))*curDT^(1) ...
                            ...+cDT2           *(DUVec.^(0))*curDT^(2) ...
                            +cDT2_noTs      *(DUVec.^(0))*curDT^(2) ...
                            +cDU3           *(DUVec.^(3))*curDT^(0) ...
                            ...+cDU2DT         *(DUVec.^(2))*curDT^(1) ...
                            +cDU2DT_noTs    *(DUVec.^(2))*curDT^(1) ...
                            ...+cDUDT2         *(DUVec.^(1))*curDT^(2) ...
                            +cDUDT2_noTs    *(DUVec.^(1))*curDT^(2) ...
                            ...+cDT3           *(DUVec.^(0))*curDT^(3) ...
                            +cDT3_noTs      *(DUVec.^(0))*curDT^(3) ...
                            +cDU4           *(DUVec.^(4))*curDT^(0) ...
                            +cDU3DT         *(DUVec.^(3))*curDT^(1) ...
                            ...+cDU2DT2        *(DUVec.^(2))*curDT^(2) ...
                            +cDU2DT2_noTs   *(DUVec.^(2))*curDT^(2) ...
                            ...+cDUDT3         *(DUVec.^(1))*curDT^(3) ...
                            +cDUDT3_noTs    *(DUVec.^(1))*curDT^(3) ...
                            ...+cDT4           *(DUVec.^(0))*curDT^(4) ...
                            +cDT4_noTs      *(DUVec.^(0))*curDT^(4) ...
                            ...+cDU6           *(DUVec.^(6))*curDT^(0) ...
                            ;
                    end
                end
            end
            figure;plot(nVec(:),rootMat);
            legend_CELL = cellfun(@(rID) ...
                [...
                '$r$=' num2str(rVec(rID)) ...
                ' limit of $\left(\frac{N}{1-r}\right)^{2}\Delta_{\theta,HPBW}$ = ' num2str(rootMat(end,rID)) ...
                ], num2cell(1:nVal_r), 'UniformOutput', false);
            legend(legend_CELL,'Interpreter','latex');
            xlabel('N');
            %ylabel('$\left(\frac{N}{1-r}\right)^{2}\Delta_{\theta,HPBW}$','Interpreter','latex');
            title('Plotting $\left(\frac{N}{1-r}\right)^{2}\Delta_{\theta,HPBW}$ vs N for various r values','Interpreter','latex');
            set(get(gca,'ylabel'),'rotation',0)
            figure;plot(rVec,rootMat(end,:));
            hold on;
            plot(rVec,2./(1-rVec));
            title({...
                'Plot of $\left(\frac{N}{1-r}\right)^{2}\Delta_{\theta,HPBW}$ vs $r$ ' ...
                ['Comparing to $\frac{2}{1-r}$'] ...
                },'Interpreter','latex');
            legend({'simulation results','$\frac{2}{1-r}$'},'Interpreter','latex');
            xlabel('$r$','Interpreter','latex');
            
            %% plot
            simResult           = real(simResult);
            [DU_MAT, DT_MAT]    = meshgrid(DUVec,DTVec);
            for n   = nVec(1:min(10,length(nVec)))
                for r   = rVec(end-5:end)
                    close all;
                    figure;
                    surf(DU_MAT, DT_MAT, simResult(:,:,nVec == n,rVec == r),    'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none');
                    hold on;
                    surf(DU_MAT, DT_MAT, approxResult(:,:,nVec == n,rVec == r),  'FaceColor','r', 'FaceAlpha',0.3, 'EdgeColor','none');
                    surf(DU_MAT, DT_MAT, approxResult_noTs(:,:,nVec == n,rVec == r),  'FaceColor','r', 'FaceAlpha',0.7, 'EdgeColor','none');
                    surf(DU_MAT, DT_MAT, 2*ones(size(DU_MAT)),                  'FaceColor','b', 'FaceAlpha',0.5, 'EdgeColor','none');
                    title(['N = ' num2str(n) ', r = ' num2str(r)]);
                    xlabel('phase [RAD]');
                    ylabel('tau [RAD]');
                    zlim([0 4]);
                    legend({'full expression' 'approx' 'approx_noTs' '3dB'});
                end
            end
        end
        
    end
    
    %% numeric validation  - with attenuation
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
end

function limVal = f_D(N,x)
limVal          = sin(N*x)./sin(x);
limVal(x==0)    = N;
end