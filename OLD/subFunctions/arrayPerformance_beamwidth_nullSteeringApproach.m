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
            syms gs positive;
            g = gs;
        end
        
        %% generating basic terms
        ad  = (g/gs)*(r/N)*exp(+1i*w*ts)*(1-exp(-1i*N*DU))/(1-exp(-1i*DU));
        bd  = ad;
        
        %% full expression transfer function
        if enable_feedback
            hDen    = 1-exp(-1i*w*t)*bd;
            h       = exp(-1i*w*t)*ad/hDen;
        else
            h       = exp(-1i*w*t)*ad;
        end
        h                   = simplify(h);
        syms dW; 
        w2  = w + dW; 
        DU2 = DU*w2/w;
        
        ad_sinCos           = simplify(rewrite(simplify(ad*exp(1i*(N-1)*DU/2)),'sincos'));
        h_sinCos            = rewrite(h,'sincos');
        ad_sideLobe         = simplify(expand(subs(ad/r,{DU r gs},{3*pi/N 0 1})));
        ad_sideLobeAbs      = simplify(rewrite(expand(ad_sideLobe*conj(ad_sideLobe)),'sincos'));
        f_ad_sideLobeAbs    = matlabFunction(ad_sideLobeAbs);
        nVec_LOCAL          = 2:100;
        %figure;plot(nVec_LOCAL,f_ad_sideLobeAbs(nVec_LOCAL));
        h_sideLobe          = subs(h,{DU r},{3*pi/N 0});
        hAbs2               = simplify(expand(h.*conj(h)));
        hAbs2_sincos        = simplify(rewrite(expand(h.*conj(h)),'sincos'));
        hAbs2_sideLobe      = simplify(expand(subs(hAbs2,{DU gs},{3*pi/N 1})));
        
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
                        try
                            hAbs2_numDiff;
                            hAbs2_denDiff;
                        catch
                            %hAbs2_numDiff = simplify(rewrite(simplify(expand(num)),'sincos'));
                            hAbs2_numDiff = simplify(num);
                            hAbs2_denDiff = simplify(den);
                        end
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
        hAbs2Rel            = simplify(rewrite(expand(hAbs2/lim_hAbs2),'sincos'));
        hAbs2_sincos_rel    = simplify(expand(hAbs2_sincos/lim_hAbs2));
        if true
            %% ambiguity
            duVecTmp            = linspace(-pi,pi,100);
            dTVecTmp            = linspace(-pi,pi,100);            
            [duMATTmp,dTMatTmp] = meshgrid(duVecTmp,dTVecTmp);
            
            hAbs2_sincos_relVal_NwdUdT = subs(hAbs2_sincos_rel,{N DU DT w},{3 duMATTmp dTMatTmp 1});
            for rTmp = [0.4 0.6 0.8 0.9]
                hAbs2_sincos_relVal = eval(subs(hAbs2_sincos_relVal_NwdUdT,{r},{rTmp}));
                figure;
                contour(duMATTmp/pi,dTMatTmp/pi,db(hAbs2_sincos_relVal),'LevelList',[0 -3 -9 -20 -40 -80],'ShowText','on');
                xticks(-0.75:0.25:0.75);
                yticks(-0.75:0.25:0.75);
            end
        end
        hAbs2Rel_onlyDU     = simplify(expand(subs(hAbs2Rel,{DT },{0})));
        if true
            %% DEBUG
            hAbs2Rel_onlyDU_values = subs(hAbs2Rel_onlyDU,{N, r, w, DU}, {5,0.7,1, linspace(-pi/2,pi/2,1000)});
            figure;plot(db(abs(hAbs2Rel_onlyDU_values)));
        end
        hAbs2Rel_DT0_gs1    = simplify(expand(rewrite(subs(hAbs2Rel,{DT gs}, {0 1}),'sincos')));
        if false
            rVec_LOCAL          = [0 0.3:0.2:0.9];
            %% sidelobe analysis
            hAbs2Rel_sideLobe   = subs(hAbs2Rel,{DU DT gs},{3*pi/N 0 1});
            hAbs2Rel_sideLobe   = simplify(expand(rewrite(hAbs2Rel_sideLobe,'sincos')));
            f_hAbs2Rel_sideLobe = matlabFunction(hAbs2Rel_sideLobe);
            nVec_LOCAL          = 2:100;
            figure;
            hold on;
            for rVal = rVec_LOCAL
                plot(nVec_LOCAL,sqrt(f_hAbs2Rel_sideLobe(nVec_LOCAL,rVal)));
                plot(nVec_LOCAL,2*(1-rVal)*ones(size(nVec_LOCAL))/(3*pi));
            end
            ylim([0,0.3]);
        end
        
        %% two freq approach
        if true
            syms DW positive;
            
            %% configure
            cVal    = 3e8;
            fS      = 1e9;
            wS      = 2*pi*fS;
            range   = 10e3;
            nVec    = 2:10;
            rVec    = [0.01 0.1:0.1:0.9 0.09];
            nValw   = 1000;
            
            %% twoFreq response
            h               = simplify(expand(h));
            w2              = w + DW;
            h2              = subs(h,w,w2);
            
            if true
                processorOutput = 1/((1/h) - (1/h2));
                twoFreqH        = simplify(rewrite(expand(processorOutput),'sincos'));
                twoFreqHAbs2    = simplify(rewrite(expand(twoFreqH*conj(twoFreqH)),'sincos'));
                if true
                    %% lim_twoFreqHAbs2_DU0
                    if true
                        curExpr         = twoFreqHAbs2;
                        [num,den]       = numden(curExpr);
                        foundExpr       = 0;
                        while ~foundExpr
                            numVal      = subs(num, DU, 0);
                            denVal      = subs(den, DU, 0);
                            if ~(denVal==0)
                                foundExpr   = 1;
                            else
                                num         = diff(num, DU);
                                den         = diff(den, DU);
                            end
                        end
                        curExpr             = numVal/denVal;
                        lim_twoFreqHAbs2_DU0= simplify(curExpr);
                    end
                    %% lim_twoFreqHAbs2_DW0
                    if true
                        curExpr         = twoFreqHAbs2;
                        [num,den]       = numden(curExpr);
                        foundExpr       = 0;
                        while ~foundExpr
                            numVal      = simplify(subs(num, DW, 0));
                            denVal      = simplify(subs(den, DW, 0));
                            if ~(denVal==0)
                                foundExpr   = 1;
                            else
                                num         = diff(num, DW);
                                den         = diff(den, DW);
                            end
                        end
                        curExpr             = numVal/denVal;
                        lim_twoFreqHAbs2_DW0= simplify(curExpr);
                    end
                end
                twoFreqHAbs2Rel     = simplify(twoFreqHAbs2 / lim_twoFreqHAbs2_DU0);
                
                phVec               = linspace(1e-6,pi,10000);
                hAbs2RelVal         = subs(hAbs2Rel,{N DT r DU}, {3 0 0.9 phVec});
                twoFreqHAbs2RelVal  = subs(twoFreqHAbs2Rel,{N DU},{3 phVec});
                figure;plot([hAbs2RelVal(:) twoFreqHAbs2RelVal(:)]);
            end
            
            %% simulate
            if false
                DUVec           = linspace(0,2*pi,100);
                DTVec           = DUVec;
                [DUMAT,DTMAT]   = meshgrid(DUVec,DTVec);
                
                parSets     = cell(length(nVec),length(rVec),length(wVec));
                for nID  = 1:length(nVec)
                    for rID=1:length(rVec)
                        for wID=1:length(wVec)
                            parSets{nID,rID,wID}.nVal       = nVec(nID);
                            parSets{nID,rID,wID}.rVal       = rVec(rID);
                            parSets{nID,rID,wID}.wVal       = wVec(wID);
                            parSets{nID,rID,wID}.DUVec      = DUVec;
                            parSets{nID,rID,wID}.DTVec      = DTVec;
                            parSets{nID,rID,wID}.DUMAT      = DUMAT;
                            parSets{nID,rID,wID}.DTMAT      = DTMAT;
                        end
                    end
                end
                simResult   = cell(length(nVec),length(rVec),length(wVec));
                for nID     = 1:length(nVec)
                    for rID     = 1:length(rVec)
                        parfor wID  = 1:length(wVec)
                            curParSet               = parSets(nID,rID,wID);
                            nVal                    = curParSet{1}.nVal;
                            rVal                    = curParSet{1}.rVal;
                            wVal                    = curParSet{1}.wVal;
                            curDUMAT                = curParSet{1}.DUMAT;
                            curDTMAT                = curParSet{1}.DTMAT;
                            curDUVec                = curParSet{1}.DUVec;
                            curDTVec                = curParSet{1}.DTVec;
                            curCRLB_DU              = subs(CRLB_DU,N,nVal);
                            curCRLB_DU              = subs(curCRLB_DU,r,rVal);
                            curCRLB_DU              = subs(curCRLB_DU,w,wVal);
                            CRLB_DU_DU0Val          = subs(CRLB_DU_DU0,DT,curDTVec);
                            curCRLB_DUnoDU0         = subs(curCRLB_DU,{DU,DT},{curDUMAT(:,2:end-1),curDTMAT(:,2:end-1)});
                            curCRLB_DUVal           = [CRLB_DU_DU0Val(:) curCRLB_DUnoDU0 CRLB_DU_DU0Val(:)];
                            simResult(nID,rID,wID)  = {curCRLB_DUVal};
                        end
                    end
                end
            end
        end
        
        %% FISHER INFORMATION MATRIX
        if false
            syms R positive;
            
            %% configure
            cVal    = 3e8;
            fS      = 1e9;
            wS      = 2*pi*fS;
            range   = 10e3;
            nVec    = 2:10;
            rVec    = [0.01 0.1:0.1:0.9 0.09];
            nValw   = 1000;
            
            %% DCBF symbolics (DT = 0)
            if true
                %% dSHAd (DTFT of ramp)
                term1   = N*exp(-1i*N*DU)*(1-exp(-1i*DU));
                term2   = (1-exp(-1i*N*DU))*exp(-1i*DU);
                term3   = 1-exp(-1i*DU);
                dSHAd   = r*exp(1i*DT)*(term1-term2)/(N*term3^2);
                %% dSHd (dirichle)
                dSHd    = r*exp(1i*DT)*(1-exp(N*DU))/(N*(1-exp(DU)));
                %% J
                denJ    = (1-dSHd)^2;
                term6   = dSHAd/denJ;
                J_DUDU  = term6*conj(term6);
                term7   = dSHd/denJ;
                J_DTDT  = w^2*term7*conj(term7);
                denJAbs = denJ*conj(denJ);
                bHdConj = dSHd;
                J_DUDT  = real(1i*w*dSHAd*bHdConj/denJAbs);
                J_DTDU  = J_DUDT;
                J       = [J_DUDU J_DUDT ; J_DTDT J_DTDU];
                %% CRLB
                CRLB    = inv(J);
                CRLB_DU = CRLB(1,1);
                CRLB_DT = CRLB(2,2);
            end
            if false
                simple = [];
                disp('real(dSHAd)');
                temp                = simplify(rewrite(real(dSHAd),'exp'))
                simple.dSHAd_real   = temp;
                disp('imag(dSHAd)');
                temp                = simplify(rewrite(imag(dSHAd),'exp'))
                simple.dSHAd_imag   = temp;
                disp('dSHAdAbs');
                temp                = simplify(rewrite(expand(dSHAd*conj(dSHAd)),'sincos'))
                simple.dSHAdAbs     = temp;
                disp('J_DUDU');
                temp                = simplify(rewrite(expand(J_DUDU),'sincos'))
                simple.J_DUDU       = temp;
                disp('J_DTDT');
                temp                = simplify(rewrite(expand(J_DTDT),'sincos'))
                simple.J_DTDT       = temp;
                disp('J_DUDT');
                temp                = simplify(rewrite(expand(J_DUDT),'sincos'))
                simple.J_DUDT       = temp;
            end
            
            %% simulate
            wVec            = linspace(-wS/2,wS/2,100);
            DUVec           = linspace(0,2*pi,100);
            DTVec           = DUVec;
            [DUMAT,DTMAT]   = meshgrid(DUVec,DTVec);
            
            if true
                %% CRLB_DU_DU0
                iterID      = 0;
                curExpr     = CRLB_DU;
                [num,den]   = numden(curExpr);
                foundExpr   = 0;
                while ~foundExpr
                    numVal      = subs(num, DU, 0);
                    denVal      = subs(den, DU, 0);
                    if ~(denVal==0)
                        foundExpr   = 1;
                    else
                        num         = diff(num, DU);
                        den         = diff(den, DU);
                        iterID      = iterID+1;
                    end
                end
                curExpr         = numVal/denVal;
                CRLB_DU_DU0     = simplify(curExpr);
                disp(['CRLB_DU_DU0 is calculated in ' num2str(iterID) ' iterations'])
                
            end
            parSets     = cell(length(nVec),length(rVec),length(wVec));
            for nID  = 1:length(nVec)
                for rID=1:length(rVec)
                    for wID=1:length(wVec)
                        parSets{nID,rID,wID}.nVal       = nVec(nID);
                        parSets{nID,rID,wID}.rVal       = rVec(rID);
                        parSets{nID,rID,wID}.wVal       = wVec(wID);
                        parSets{nID,rID,wID}.DUVec      = DUVec;
                        parSets{nID,rID,wID}.DTVec      = DTVec;
                        parSets{nID,rID,wID}.DUMAT      = DUMAT;
                        parSets{nID,rID,wID}.DTMAT      = DTMAT;
                    end
                end
            end
            simResult   = cell(length(nVec),length(rVec),length(wVec));
            for nID     = 1:length(nVec)
                for rID     = 1:length(rVec)
                    parfor wID  = 1:length(wVec)
                        curParSet               = parSets(nID,rID,wID);
                        nVal                    = curParSet{1}.nVal;
                        rVal                    = curParSet{1}.rVal;
                        wVal                    = curParSet{1}.wVal;
                        curDUMAT                = curParSet{1}.DUMAT;
                        curDTMAT                = curParSet{1}.DTMAT;
                        curDUVec                = curParSet{1}.DUVec;
                        curDTVec                = curParSet{1}.DTVec;
                        curCRLB_DU              = subs(CRLB_DU,N,nVal);
                        curCRLB_DU              = subs(curCRLB_DU,r,rVal);
                        curCRLB_DU              = subs(curCRLB_DU,w,wVal);
                        CRLB_DU_DU0Val          = subs(CRLB_DU_DU0,DT,curDTVec);
                        curCRLB_DUnoDU0         = subs(curCRLB_DU,{DU,DT},{curDUMAT(:,2:end-1),curDTMAT(:,2:end-1)});
                        curCRLB_DUVal           = [CRLB_DU_DU0Val(:) curCRLB_DUnoDU0 CRLB_DU_DU0Val(:)];
                        simResult(nID,rID,wID)  = {curCRLB_DUVal};
                    end
                end
            end
        end
        
        %% taylor
        if false
            %% d/dDU : 0
            hAbs2Rel_dDU        = diff(hAbs2Rel,DU);
            if true
                curExpr     = simplify(subs(hAbs2Rel_dDU,DT,0));
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
            
            %% d/dDU2 : -\left(\left(N - 1\right)*\left(N - 4*r + 2*N*r + 1\right)\right)/\left(6*\left(r - 1\right)^2\right)
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
            
            %% d/dDUDT : -\left(r*w*\left(N - 1\right)\right)/\left(r - 1\right)^2
            hAbs2Rel_dDUDT  = diff(hAbs2Rel_dDT,DU);
            if true
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
            
            %% d/dDT2 : -\left(2*r*w^2\right)/\left(r - 1\right)^2
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
            if true
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
            if true
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
            if true
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
            if true
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
            
            %% d/dDU4 : -\left(\left(N - 1\right)*\left(2*N^3*r^3 - 21*N^3*r^2 - 24*N^3*r - 2*N^3 + 2*N^2*r^3 + 99*N^2*r^2 + 36*N^2*r - 2*N^2 - 18*N*r^3 - 126*N*r^2 + 6*N*r + 3*N + 12*r^3 + 54*r^2 - 24*r + 3\right)\right)/\left(30*\left(r - 1\right)^4\right)
            hAbs2Rel_dDU4   = diff(hAbs2Rel_dDU3,DU);
            if true
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
            
            %% d/dDT3DU : \left(r*w^3*\left(N - 1\right)*\left(r^2 + 10*r + 1\right)\right)/\left(r - 1\right)^4
            hAbs2Rel_dDT3DU = diff(hAbs2Rel_dDT3,DU);
            if true
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
            
            %% d/dDU2DT2 : \left(r*w^2*\left(N - 1\right)*\left(3*N - 16*r + 14*N*r + N*r^2 - 2*r^2\right)\right)/\left(3*\left(r - 1\right)^4\right)
            hAbs2Rel_dDU2dDT2  = diff(hAbs2Rel_dDUDT2,DU);
            if true
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
            
            %% d/dDTDU3 : \left(r*w*\left(N - 1\right)^2*\left(2*N - 6*r + 4*N*r - r^2 + 1\right)\right)/\left(2*\left(r - 1\right)^4\right)
            hAbs2Rel_dDTDU3 = diff(hAbs2Rel_dDU3,DT);
            if true
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
            
            %% d/dDT4 : \left(2*r*w^4*\left(r^2 + 10*r + 1\right)\right)/\left(r - 1\right)^4
            hAbs2Rel_dDT4   = diff(hAbs2Rel_dDT3,DT);
            if true
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
            if true
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
            if true
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
            if true
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
            if true
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
            if true
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
            if true
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
                init_N      = 5;
                final_N     = 100;
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
            rVec    = [0 0.1 : 0.05 : 0.95 , 0.99];%DO NOT CHANGE old = linspace(init_r,final_r,nVal_r);
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
            rootMat_fminbnd     = zeros(nVal_N,length(rVec));
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
            uVec            = linspace(1e-6,pi,100);
            
            for n = nVec
                symVarValues(2) = n;
                for r = rVec
                    symVarValues(4) = r;
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
                    
                    f_curBp             = @(pVec) (f_hAbs2Rel(symVarValues,pVec(:)')-0.5).^2;
                    
                    distanceVec         = f_curBp(uVec);
                    [~,minDistanceID]   = min(distanceVec);
                    minDistanceU        = uVec(minDistanceID);
                    rootFMinBnd         = fminbnd(f_curBp,0.9*minDistanceU,1.1*minDistanceU);
                    curRoot             = rootFMinBnd;
                    fEvalAtRoot         = f_hAbs2Rel(symVarValues,curRoot);
                    if true
                        close all;
                        phasValues          = logspace(-6,log10(pi),100000);
                        symVarValuesULA     = symVarValues;
                        symVarValuesULA(4)  = 0;
                        ulaBP               = f_hAbs2Rel(symVarValuesULA,phasValues);
                        ulaBP(isnan(ulaBP)) = 1;
                        beamPattern         = f_hAbs2Rel(symVarValues,phasValues);
                        beamPattern(isnan(beamPattern)) = 1;
                        beamPatternAux      = f_curBp(phasValues);
                        beamPatternApprox2  = 1 + cDU2*phasValues.^2;
                        beamPatternApprox4  = 1 + cDU2*phasValues.^2 + cDU4*phasValues.^4;
                        beamPatternApprox6  = 1 + cDU2*phasValues.^2 + cDU4*phasValues.^4 + cDU6*phasValues.^6;
                        figure;
                        %                         plot(phasValues,...
                        %                             db([...
                        %                             ulaBP(:)...
                        %                             ...beamPattern(:) ...
                        %                             ...beamPatternAux(:) ...
                        %                             ...beamPatternApprox2(:) ...
                        %                             ...beamPatternApprox4(:) ...
                        %                             ...beamPatternApprox6(:)...
                        %                             ]) ...
                        %                             ,'*'...
                        %                             );
                        %                         hold on;
                        plot(phasValues,...
                            db([...
                            ...ulaBP(:)...
                            beamPattern(:) ...
                            ...beamPatternAux(:) ...
                            beamPatternApprox2(:) ...
                            beamPatternApprox4(:) ...
                            beamPatternApprox6(:)...
                            ]));
                        legend('Beampattern','Approx2','Approx4','Approx6');
                        title(['N=' num2str(n) ', r=' num2str(r)]);
                        ylim([-4,.5]);
                    end
                    rootMat(...
                        nVec    == n,       ... n
                        rVec    == r        ... r
                        ) = curRoot6*n/2;
                    rootMat_fminbnd(...
                        nVec    == n,       ... n
                        rVec    == r        ... r
                        ) = rootFMinBnd*n/2;
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
            plotRIDVec  = [1 2 : 4 : length(rVec)-1];
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
            polyCoeffs      = polyfit(plotRVec,rootMat(end,plotRIDVec)./(plotRVec-1),2);
            syms x;
            polyCoeffsRound = round(10*polyCoeffs)/10;
            residualPoly    = (polyCoeffsRound*reshape(x.^(2:-1:0),[],1));
            polySym         = (x - 1).*residualPoly;
            polyValue       = eval(subs(polySym,x,plotRVec));
            plot(plotRVec,polyValue,'o-','MarkerIndices',1:2:length(plotRVec));
            ylim([0 1.5])
            %             title({...
            %                 'Plot of $\right))(\frac{N}{1-r}\right)^{2}\Delta_{\theta,HPBW}$ vs $r$ ' ...
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