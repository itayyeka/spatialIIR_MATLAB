function [simBp_dbAbs_norm]     = calc_simBp_dbAbs_norm(simOut)
simBp_dbAbs             = db(abs(simOut.hMat(end,:)).^2);
simBp_dbAbs_norm        = simBp_dbAbs-max(simBp_dbAbs(:));
end