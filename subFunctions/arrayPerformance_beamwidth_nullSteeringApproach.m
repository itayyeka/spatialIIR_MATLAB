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
        
        f_hAbs2 = matlabFunction(hAbs2,'Vars',{ [DT, N, c, r, ts, w]' DU});
        
        if true
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
            
            %% d/du
            hAbs2Rel        = lim_hAbs2/hAbs2;
            hAbs2Rel_dDU    = diff(hAbs2Rel,DU);
            if true
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
            
            %% d/du2
            hAbs2Rel_dDU2   = diff(hAbs2Rel_dDU,DU);
            if true
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
            
            %% d/dDT
            hAbs2Rel_dDT    = diff(hAbs2Rel,DT);
            if true
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
            
            %% d/dDT2
            hAbs2Rel_dDT2 = diff(hAbs2Rel_dDT,DT);
            if true
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
            
            %% d/dudDT            
            hAbs2Rel_dDUDT  = diff(hAbs2Rel_dDT,DU);
            if true
                syms R phi;
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
        end
        
        %% numeric validation
        if false
            %% configure
            if true
                %% N
                init_N      = 4;
                final_N     = 10;
                %% r
                nVal_r      = 10;
                init_r      = 0.6;
                final_r     = 0.99;
                %% u
                nVal_u      = 10;
                init_u      = 1e-6;
                final_u     = 0.2;
                %% t
                nVal_t      = 10;
                init_t      = 0;
                final_t     = 0.2;
                %% DF
                fVal        = 10e9;
                nVal_DF     = 10;
                init_DF     = fVal/1000;
                final_DF    = fVal/1000000;
            end
            
            %% generate input
            nVec    = init_N : final_N;
            rVec    = linspace(init_r,final_r,nVal_r);
            uVec    = linspace(init_u,final_u,nVal_u);
            tVec    = linspace(init_t,final_t,nVal_t);
            DFVec   = linspace(init_DF,final_DF,nVal_DF);
            
            %% simulate
            nVal_N          = length(nVec);
            simResult       = zeros(nVal_DF,nVal_t,nVal_u,nVal_N,nVal_r);
            symVarValues    = zeros(5,1);
            
            symVarValues(5) = fVal;
            for r = rVec
                symVarValues(1) = r;
                for n = nVec
                    symVarValues(2) = n;
                    for t = tVec
                        symVarValues(3) = t;
                        for DF = DFVec
                            symVarValues(4) = DF;
                            simResult(...
                                DFVec   == DF,  ... DF
                                tVec    == t,   ... t
                                :,              ... u
                                nVec    == n,   ... n
                                rVec    == r    ... r
                                ) = ...
                                f_hOrig(symVarValues,uVec);
                        end
                    end
                end
            end
            symResultAbs2 = simResult.*conj(simResult);
            
            %% plot
            if true
                %% responseMag vs DF & tau
                plotData    = squeeze( ...
                    symResultAbs2(...
                    :,      ... DF
                    :,      ... t
                    1,      ... u
                    1,      ... n
                    end     ... r
                    ));
                
                [DF_MAT, t_MAT]     = meshgrid(DFVec,tVec);
                figure; surf(DF_MAT,t_MAT,db(plotData));
            end
        end
        
        %% Taylor expansion
        if false
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
        end
        
    end
    
    %% numeric validation  - with attenuation
end

function limVal = f_D(N,x)
limVal          = sin(N*x)./sin(x);
limVal(x==0)    = N;
end