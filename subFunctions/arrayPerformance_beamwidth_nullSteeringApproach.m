clear all;
close all;
clc;


if true
    %% symbolic calc - with attenuation
    if true
        %% configure
        enable_feedback                     = 1;
        enable_attenuationSpatialDependency = 0;
        
        %% symbolics
        syms N DU ts DT c w r positive;
        syms k;
        
        t   = ts + DT;
        
        if enable_attenuationSpatialDependency
            gs  = 1/(c*ts);
            g   = 1/(c*t);
        else
            syms gs;
            g = gs;
        end
        
        %% generating basic terms
        if enable_feedback
            ad  = (r/N)*(1-exp(-1i*N*DU))./(1-exp(-1i*DU));%exp(-1i*(DU*(N-1)/2)).*sin(N*DU/2)./sin(DU/2);
        else
            ad  = (1/N)*(1-exp(-1i*N*DU))./(1-exp(-1i*DU));%exp(-1i*(DU*(N-1)/2)).*sin(N*DU/2)./sin(DU/2);
        end
        bd  = ad * exp(1i*w*ts) / (gs^2) ;
        
        %% full expression transfer function
        if enable_feedback
            h   = (ad*(gs^2))./(1-(g^2)*bd*exp(-1i*w*t));
        else
            h   = ad;
        end
        hAbs2   = h.*conj(h);
        
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
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                catch
                    num         = diff(num, DU);
                    den         = diff(den, DU);
                end
            end
            curExpr     = numVal/denVal;
            lim_hAbs2   = simplify(curExpr)
        end
        hAbs2Rel        = simplify(expand(rewrite(hAbs2/lim_hAbs2,'sincos')));
        
        %% taylor
        if false
            %% d/dDU : 0
            hAbs2Rel_dDU   = diff(hAbs2Rel,DU);
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
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU   = simplify(curExpr)
            end
            
            %% d/dDT : 0
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
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDT   = simplify(curExpr)
            end
            
            %% d/dDU2 : -((N - 1)*(N - 4*r + 2*N*r + 1))/(6*(r - 1)^2)
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
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU2   = simplify(curExpr)
            end
            
            %% d/dDUDT : -(r*w*(N - 1))/(r - 1)^2
            hAbs2Rel_dDUDT  = diff(hAbs2Rel_dDT,DU);
            if false
                curExpr     = subs(hAbs2Rel_dDUDT,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDUDT  = simplify(curExpr)
            end
            
            %% d/dDT2 : -(2*r*w^2)/(r - 1)^2
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
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
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
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU3   = simplify(curExpr)
            end
            
            %% d/dDU2DT : 0
            hAbs2Rel_dDU2DT = diff(hAbs2Rel_dDU2,DT);
            if false
                curExpr     = subs(hAbs2Rel_dDU2DT,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU2DT  = simplify(curExpr)
            end
            
            %% d/dDUDT2 : 0
            hAbs2Rel_dDUDT2 = diff(hAbs2Rel_dDT2,DU);
            if false
                curExpr     = subs(hAbs2Rel_dDUDT2,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDUDT2  = simplify(curExpr)
            end
            
            %% d/dDT3 : 0
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
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDT3   = simplify(expand(curExpr))
            end
            
            %% d/dDU4 : -((N - 1)*(2*N^3*r^3 - 21*N^3*r^2 - 24*N^3*r - 2*N^3 + 2*N^2*r^3 + 99*N^2*r^2 + 36*N^2*r - 2*N^2 - 18*N*r^3 - 126*N*r^2 + 6*N*r + 3*N + 12*r^3 + 54*r^2 - 24*r + 3))/(30*(r - 1)^4)
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
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU4   = simplify(curExpr)
            end
            
            %% d/dDT3DU : (r*w^3*(N - 1)*(r^2 + 10*r + 1))/(r - 1)^4
            hAbs2Rel_dDT3DU = diff(hAbs2Rel_dDT3,DU);
            if false
                curExpr     = subs(hAbs2Rel_dDT3DU,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDT3DU  = simplify(curExpr)
            end
            
            %% d/dDU2DT2 : (r*w^2*(N - 1)*(3*N - 16*r + 14*N*r + N*r^2 - 2*r^2))/(3*(r - 1)^4)
            hAbs2Rel_dDU2dDT2  = diff(hAbs2Rel_dDUDT2,DU);
            if false
                curExpr     = subs(hAbs2Rel_dDU2dDT2,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU2dDT2  = simplify(curExpr)
            end
            
            %% d/dDTDU3 : (r*w*(N - 1)^2*(2*N - 6*r + 4*N*r - r^2 + 1))/(2*(r - 1)^4)
            hAbs2Rel_dDTDU3 = diff(hAbs2Rel_dDU3,DT);
            if false
                curExpr     = subs(hAbs2Rel_dDTDU3,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDTDU3  = simplify(curExpr)
            end
            
            %% d/dDT4 : (2*r*w^4*(r^2 + 10*r + 1))/(r - 1)^4
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
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
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
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU5   = simplify(curExpr)
            end
            
            %% d/dDU4DT : 0
            hAbs2Rel_dDU4DT = diff(hAbs2Rel_dDU4,DT);
            if false
                curExpr     = subs(hAbs2Rel_dDU4DT,DT, 0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU4DT  = simplify(curExpr)
            end
            
            %% d/dDU3DT2 : 0
            hAbs2Rel_dDU3DT2 = diff(hAbs2Rel_dDU2dDT2,DU);
            if false
                curExpr     = subs(hAbs2Rel_dDU3DT2,DT, 0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = diff(num, DU);
                            den         = (diff(den, DU));
                        end
                    catch
                        num         = diff(num, DU);
                        den         = (diff(den, DU));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU3DT2  = simplify(curExpr)
            end
            
            %% d/dDU2DT3 : 0
            hAbs2Rel_dDU2DT3 = diff(hAbs2Rel_dDU2dDT2,DT);
            if false
                curExpr     = subs(hAbs2Rel_dDU2DT3,DT, 0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = diff(num, DU);
                            den         = (diff(den, DU));
                        end
                    catch
                        num         = diff(num, DU);
                        den         = (diff(den, DU));
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU2DT3  = simplify(curExpr)
            end
            
            %% d/dDUDT4 : 0
            hAbs2Rel_dDUDT4 = diff(hAbs2Rel_dDT3DU,DT);
            if false
                curExpr     = subs(hAbs2Rel_dDUDT4,DT, 0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDUDT4  = simplify(curExpr)
            end
            
            %% d/dDT5 : 0
            hAbs2Rel_dDT5   = diff(hAbs2Rel_dDT4,DT);
            if false
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
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDT5   = simplify(expand(curExpr))
            end
            
            %% d/dDU6 : -((N - 1)*(4*N^5*r^5 - 62*N^5*r^4 + 264*N^5*r^3 + 590*N^5*r^2 + 146*N^5*r + 3*N^5 + 4*N^4*r^5 + 106*N^4*r^4 - 2340*N^4*r^3 - 2266*N^4*r^2 - 232*N^4*r + 3*N^4 - 24*N^3*r^5 + 645*N^3*r^4 + 6200*N^3*r^3 + 2760*N^3*r^2 - 120*N^3*r - 11*N^3 - 24*N^2*r^5 - 1875*N^2*r^4 - 7240*N^2*r^3 - 600*N^2*r^2 + 300*N^2*r - 11*N^2 + 60*N*r^5 + 1779*N*r^4 + 3848*N*r^3 - 978*N*r^2 + 6*N*r + 10*N - 24*r^5 - 573*r^4 - 772*r^3 + 534*r^2 - 120*r + 10))/(84*(r - 1)^6)
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
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU6   = simplify(curExpr)
            end
            
            %% d/dDU7 : 
            hAbs2Rel_dDU7   = diff(hAbs2Rel_dDU6,DU);
            if true
                curExpr     = subs(hAbs2Rel_dDU7,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU7   = simplify(curExpr)
            end
            
            %% d/dDU8 : 
            hAbs2Rel_dDU8   = diff(hAbs2Rel_dDU7,DU);
            if true
                curExpr     = subs(hAbs2Rel_dDU8,DT,0);
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, DU, 0);
                        denVal      = subs(den, DU, 0);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = diff(num, DU);
                            den         = diff(den, DU);
                        end
                    catch
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                end
                curExpr             = numVal/denVal;
                lim_hAbs2Rel_dDU8   = simplify(curExpr)
            end
            
        end
        
        %% numeric validation
        if enable_attenuationSpatialDependency
            symVarVec = [DT, N, c, r, ts, w]';
        else
            symVarVec = transpose([DT, N, c, r, ts, w, gs]);
        end
        
        f_hAbs2Rel  = matlabFunction(hAbs2Rel,'Vars',{symVarVec DU});
        if true
            %% configure
            cVal        = 3e8;
            fVal        = 10e9;
            RVal        = 10e3;
            lambdaVal   = cVal/fVal;
            if true
                %% N
                init_N      = 2;
                final_N     = 200;
                hop_N       = 1;
                %% r
                nVal_r      = 30;
                init_r      = 0.;
                final_r     = 0.99;
                %% DU
                nVal_DU     = 40;
                init_DU     = 1e-12;
                final_DU    = 1;
                %% DT
                nVal_DT     = 30;
                init_DT     = 0;
                final_DT    = .04*lambdaVal/cVal;
            end            
            %% generate input
            nVec    = init_N : hop_N : final_N;
            rVec    = [0 0.1 : 0.05 : 0.95 , 0.99];%linspace(init_r,final_r,nVal_r);
            DUVec   = linspace(init_DU,final_DU,nVal_DU);
            DTVec   = linspace(init_DT,final_DT,nVal_DT);            
            if ~enable_feedback
                rVec = 0;
            end            
            %% simulate
            nVal_N              = length(nVec);
            simResult           = zeros(nVal_DT,nVal_DU,nVal_N,length(rVec));
            approxResult        = zeros(nVal_DT,nVal_DU,nVal_N,length(rVec));
            approxResult_noTs   = zeros(nVal_DT,nVal_DU,nVal_N,length(rVec));
            rootMat             = zeros(nVal_N,length(rVec));
            symVarValues        = zeros(7,1);
            
            %[DT, N, c, r, ts, w, gs]
            wVal    = 2*pi*fVal;
            tsVal   = RVal/cVal;
            if ~enable_attenuationSpatialDependency
                gsVal   = 1;
                symVarValues(7) = gsVal;
            end
            symVarValues(6) = wVal;
            symVarValues(3) = cVal;
            symVarValues(5) = tsVal;
            
            for r = rVec
                symVarValues(4) = r;
                for n = nVec
                    symVarValues(2) = n;
                    cDU             = (nchoosek(1,0)/factorial(1))*(0);
                    cDT             = (nchoosek(1,1)/factorial(1))*(0);
                    cDU2            = (nchoosek(2,0)/factorial(2))*(-((n - 1)*(n - 4*r + 2*n*r + 1))/(6*(r - 1)^2));
                    cDUDT           = (nchoosek(2,1)/factorial(2))*(-(r*wVal*(n - 1))/(r - 1)^2);
                    cDT2            = (nchoosek(2,2)/factorial(2))*(-(2*r*wVal^2)/(r - 1)^2);
                    cDT2_noTs       = (nchoosek(2,2)/factorial(2))*(0);
                    cDU3            = (nchoosek(3,0)/factorial(3))*(0);
                    cDU2DT          = (nchoosek(3,1)/factorial(3))*(0);
                    cDUDT2          = (nchoosek(3,2)/factorial(3))*(0);
                    cDT3            = (nchoosek(3,3)/factorial(3))*(0);
                    cDU4            = (nchoosek(4,0)/factorial(4))*(-((n - 1)*(2*n^3*r^3 - 21*n^3*r^2 - 24*n^3*r - 2*n^3 + 2*n^2*r^3 + 99*n^2*r^2 + 36*n^2*r - 2*n^2 - 18*n*r^3 - 126*n*r^2 + 6*n*r + 3*n + 12*r^3 + 54*r^2 - 24*r + 3))/(30*(r - 1)^4));
                    cDU3DT          = (nchoosek(4,1)/factorial(4))*((r*wVal*(n - 1)^2*(2*n - 6*r + 4*n*r - r^2 + 1))/(2*(r - 1)^4));
                    cDU2DT2         = (nchoosek(4,2)/factorial(4))*((r*wVal^2*(n - 1)*(3*n - 16*r + 14*n*r + n*r^2 - 2*r^2))/(3*(r - 1)^4));
                    cDUDT3          = (nchoosek(4,3)/factorial(4))*((r*wVal^3*(n - 1)*(r^2 + 10*r + 1))/(r - 1)^4);
                    cDT4            = (nchoosek(4,4)/factorial(4))*((2*r*wVal^4*(r^2 + 10*r + 1))/(r - 1)^4);
                    cDU5            = (nchoosek(5,0)/factorial(5))*(0);
                    cDU4DT          = (nchoosek(5,1)/factorial(5))*(0);
                    cDU3DT2         = (nchoosek(5,2)/factorial(5))*(0);
                    cDU2DT3         = (nchoosek(5,3)/factorial(5))*(0);
                    cDUDT4          = (nchoosek(5,4)/factorial(5))*(0);
                    cDT5            = (nchoosek(5,5)/factorial(5))*(0);
                    cDU6            = (nchoosek(6,0)/factorial(6))*(-((n - 1)*(4*n^5*r^5 - 62*n^5*r^4 + 264*n^5*r^3 + 590*n^5*r^2 + 146*n^5*r + 3*n^5 + 4*n^4*r^5 + 106*n^4*r^4 - 2340*n^4*r^3 - 2266*n^4*r^2 - 232*n^4*r + 3*n^4 - 24*n^3*r^5 + 645*n^3*r^4 + 6200*n^3*r^3 + 2760*n^3*r^2 - 120*n^3*r - 11*n^3 - 24*n^2*r^5 - 1875*n^2*r^4 - 7240*n^2*r^3 - 600*n^2*r^2 + 300*n^2*r - 11*n^2 + 60*n*r^5 + 1779*n*r^4 + 3848*n*r^3 - 978*n*r^2 + 6*n*r + 10*n - 24*r^5 - 573*r^4 - 772*r^3 + 534*r^2 - 120*r + 10))/(84*(r - 1)^6));
                    cDU8            = (nchoosek(8,0)/factorial(8))*(0);
                    cVec6           = [cDU6 cDU4 cDU2 1-0.5]; % compare to 0.5
                    rootVec6        = roots(cVec6);
                    curRoot6        = sqrt(rootVec6(cellfun(@isreal,(num2cell(rootVec6)))));
                    curRoot         = curRoot6;
                    fEvalAtRoot     = f_hAbs2Rel(symVarValues,curRoot);
                    if false
                        phasValues          = linspace(1e-6,2,100);
                        beamPattern         = f_hAbs2Rel(symVarValues,phasValues);
                        beamPatternApprox2  = 1 + cDU2*phasValues.^2;
                        beamPatternApprox4  = 1 + cDU2*phasValues.^2 + cDU4*phasValues.^4;
                        beamPatternApprox6  = 1 + cDU2*phasValues.^2 + cDU4*phasValues.^4 + cDU6*phasValues.^6;
                        figure;plot(phasValues,[beamPattern(:) beamPatternApprox2(:) beamPatternApprox4(:) beamPatternApprox6(:)]);
                        legend('Beampattern','Approx2','Approx4','Approx6');
                        ylim([0,1.5]);
                    end
                    rootMat(...
                        nVec    == n,       ... n
                        rVec    == r        ... r
                        ) = curRoot*n/2;
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
                            ...+cDU3           *(DUVec.^(3))*curDT^(0) ...
                            ...+cDU2DT         *(DUVec.^(2))*curDT^(1) ...
                            ...+cDU2DT_noTs    *(DUVec.^(2))*curDT^(1) ...
                            ...+cDUDT2         *(DUVec.^(1))*curDT^(2) ...
                            ...+cDUDT2_noTs    *(DUVec.^(1))*curDT^(2) ...
                            ...+cDT3           *(DUVec.^(0))*curDT^(3) ...
                            ...+cDT3_noTs      *(DUVec.^(0))*curDT^(3) ...
                            ...+cDU4           *(DUVec.^(4))*curDT^(0) ...
                            ...+cDU3DT         *(DUVec.^(3))*curDT^(1) ...
                            ...+cDU2DT2        *(DUVec.^(2))*curDT^(2) ...
                            ...+cDU2DT2_noTs   *(DUVec.^(2))*curDT^(2) ...
                            ...+cDUDT3         *(DUVec.^(1))*curDT^(3) ...
                            ...+cDUDT3_noTs    *(DUVec.^(1))*curDT^(3) ...
                            ...+cDT4           *(DUVec.^(0))*curDT^(4) ...
                            ...+cDT4_noTs      *(DUVec.^(0))*curDT^(4) ...
                            ...+cDU6           *(DUVec.^(6))*curDT^(0) ...
                            ;                        
                    end
                end
            end
            plotNIDVec  = 1 : length(nVec);
            plotRIDVec  = [1 2 : 5 : length(rVec)-1, length(rVec)];
            plotNVec    = nVec(plotNIDVec);
            plotRVec    = rVec(plotRIDVec);
            figure;plot(nVec,rootMat(:,plotRIDVec));
            AAA=1;
            legend_CELL = cellfun(@(rID) ...
                [...
                '$r$=' num2str(plotRVec(rID)) ...
                '\ $\lim{}N\Delta_{\theta,HPBW}$ = ' num2str(rootMat(end,plotRIDVec(rID))) ...
                ], num2cell(1:length(plotRVec)), 'UniformOutput', false);
            legend(legend_CELL,'Interpreter','latex');
            xlabel('N','Interpreter','latex');
            ylabel('$N\Delta_{\theta,HPBW}$','Interpreter','latex');
            %title('Plotting $N\Delta_{\theta,HPBW}$ vs N for various r values','Interpreter','latex');
            set(get(gca,'ylabel'),'rotation',0)
            %rootMat = zeros(nVal_N,nVal_r);
            plotNIDVec  = 1 : 5 : length(nVec);
            plotRIDVec  = 1 : length(rVec);
            plotNVec    = nVec(plotNIDVec);
            plotRVec    = rVec(plotRIDVec);
            figure;plot(plotRVec,rootMat(end,plotRIDVec),'square-','MarkerIndices',1:2:length(plotRVec));
            hold on;
            polyValues  = (rVec.^2)-2.4*(rVec.^1)+1.4;
            plot(rVec,polyValues,'O-','MarkerIndices',1:2:length(plotRVec));
            
            %             title({...
            %                 'Plot of $\left(\frac{N}{1-r}\right)^{2}\Delta_{\theta,HPBW}$ vs $r$ ' ...
            %                 ['Comparing to $\frac{2}{1-r}$'] ...
            %                 },'Interpreter','latex');
            [hleg1, hobj1] = legend({'' ''},'Interpreter','latex','position',[0.45 0.75 0.2 0.15]);
            legend boxoff;
            %             xlabel('$r$','Interpreter','latex');
            %             ylabel('$\lim{}N\Delta_{\theta,HPBW}$','Interpreter','latex');
            %% plot
            simResult           = real(simResult);
            [DU_MAT, DT_MAT]    = meshgrid(DUVec,DTVec);
            if true
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
                        %simResult           = zeros(nVal_DT,nVal_DU,nVal_N,nVal_r);
                        %                         figure;plot(DUVec,simResult(1,:,nVec == n,rVec == r));
                        %                         ylim([0 4]);
                    end
                end
            end
            if true
                plotNIDVec = 1 : 5 : length(nVec);
                plotRIDVec = 1 : 5 : length(rVec);
                plotNVec    = nVec(plotNIDVec);
                plotRVec    = rVec(plotRIDVec);
                legendCELL  = cellfun(@(rID) ['$r = ' num2str(plotRVec(rID)) '$'], num2cell(1:length(plotRVec)), 'UniformOutput', false);
                for n   = plotNVec
                    close all;
                    figure;
                    hold on;
                    for r   = plotRVec
                        plot(DUVec,simResult(1,:,nVec == n,rVec == r));
                    end
                    ylim([0.9 3]);
                    plot(DUVec,2*ones(size(DUVec)));
                    legend(legendCELL,'Interpreter','latex');
                end
            end
        end
    end
    
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
        syms DU w positive;
        hOrigAbs_d_dudt  = simplify(subs(diff(hOrigAbs_d_dt,u),{u t},{DU*cos(w) DU*sin(w)}));
        if false
            curExpr     = hOrigAbs_d_dudt;
            foundExpr   = 0;
            [num,den]   = numden(curExpr);
            while ~foundExpr
                try
                    numVal  = subs(num, DU, 0);
                    denVal  = subs(den, DU, 0);
                    if ~(denVal==0)
                        foundExpr   = 1;
                    else
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                    end
                catch
                    num         = diff(num, DU);
                    den         = diff(den, DU);
                end
            end
            curExpr             = simplify(numVal/denVal);
            lim_hOrigAbs_d_dudt = curExpr
        end
    end
end

function limVal = f_D(N,x)
limVal          = sin(N*x)./sin(x);
limVal(x==0)    = N;
end