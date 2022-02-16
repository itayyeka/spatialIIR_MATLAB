function [simOut] = spatialIIR_dualFreq(simCfg)
    standaloneFalg  = 0;
    try
        simCfg;
    catch
        close all;
        standaloneFalg = 1;
        simCfg = f_get_dualFreqCfg();
    end
    simCfg_base     = simCfg;
    N               = simCfg.nSensors;

    if true
        %% ULA
        if standaloneFalg
            simCfg          = simCfg_base;
            simCfg.r        = 0;
            simOut_ULA      = spatialIIR_singleFreq(simCfg);
            bp_ULA          = db(abs(simOut_ULA.hMat(end,:)).^2);
            bp_ULA_norm     = bp_ULA-max(bp_ULA(:));
        end
        %% singleFreq
        if standaloneFalg
            simCfg              = simCfg_base;
            simOut_singleFreq   = spatialIIR_singleFreq(simCfg);
            bp_singleFreq       = db(abs(simOut_singleFreq .hMat(end,:)).^2);
            bp_singleFreq_norm  = bp_singleFreq-max(bp_singleFreq(:));
        end
        %% f1
        if true
            simCfg                      = simCfg_base;
            simCfg.r                    = simCfg_base.r*simCfg_base.kappa;
            simCfg.compensationFreq     = 0;%-simCfg_base.sigFreq*simCfg.dF_ratio;
            simOut_f1                   = spatialIIR_singleFreq(simCfg);
            stft1                       = simOut_f1.stftVec;
        end
        %% f2
        if true
            simCfg                          = simCfg_base;
            simCfg.overrideAlpha            = [1 zeros(1,N-1)];
            simCfg.overrideBeta             = -simCfg.overrideAlpha;
            simCfg.sigFreq                  = simCfg_base.sigFreq*(1+simCfg.dF_ratio);
            simOut_f2                       = spatialIIR_singleFreq(simCfg);
            stft2                           = simOut_f2.stftVec;
            hTwoFreq_err0                   = 1./(1./stft1 + 1./stft2);
            hTwoFreq_err0_dbAbs2            = db(abs(hTwoFreq_err0).^2);
            hTwoFreq_err0_dbAbs2_norm       = hTwoFreq_err0_dbAbs2 - max(hTwoFreq_err0_dbAbs2(:));
        end
        %% theory
        targetAngleVec          = simOut_f1.cfg.targetAngleVec;
        theoryBp_dbAbs_norm     = f_theoryBp_dbAbs_norm(simOut_f1.cfg.targetDuVec,simCfg);
        nTheta                  = simOut_f1.cfg.nTheta;
        %% simOut
        if standaloneFalg
            simOut.theoryBp_dbAbs_norm  = theoryBp_dbAbs_norm;
            simOut.bp_singleFreq_norm   = bp_singleFreq_norm;
        end
        simOut.hTwoFreq_err0_dbAbs2_norm    = hTwoFreq_err0_dbAbs2_norm;
        simOut.targetAngleVec               = targetAngleVec;
        simOut.nTheta                       = nTheta;
        simOut.cfg                          = simOut_f1.cfg;

        %% plot
        if standaloneFalg
            figure;
            hold on;
            plot(targetAngleVec,theoryBp_dbAbs_norm(:),'*b-','MarkerIndices',1:6:nTheta);
            plot(targetAngleVec,bp_ULA_norm,'.-','MarkerIndices',3:6:nTheta);
            plot(targetAngleVec,bp_singleFreq_norm,'squarer-','MarkerIndices',3:6:nTheta);
            plot(targetAngleVec,hTwoFreq_err0_dbAbs2_norm(:),'diamondg-','MarkerIndices',5:6:nTheta);
        end
    end
end