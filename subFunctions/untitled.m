%% symbolic calc
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
                curExpr                 = numVal/denVal;
                lim_hAbs2Rel_du   = simplify(curExpr)
            end
            
            %% d/du2
            hAbs2Rel_dDU2  = diff(hAbs2Rel_dDU,DU);
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
                curExpr                 = numVal/denVal;
                lim_hAbs2Rel_dDU2  = simplify(curExpr)
            end
            
            %% d/dDT
            hAbs2Rel_dDT  = diff(hAbs2Rel,DT);
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
                curExpr                 = numVal/denVal;
                lim_hAbs2Rel_dDT   = simplify(curExpr)
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
                curExpr                 = numVal/denVal;
                lim_hAbs2Rel_dDUDT = simplify(curExpr)
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