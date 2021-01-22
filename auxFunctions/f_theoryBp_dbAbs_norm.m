function [theoryBp_dbAbs_norm]  = f_theoryBp_dbAbs_norm(duVec,r,N)
f_theoryBp          = @(x,r,N) ...
    (r^2)*(cos(N*x)-1) ...
    ./ ...
    (...
    N^2*(cos(x)-1) ...
    +r^2*(cos(N*x)-1) ...
    +N*r*(1+cos((N-1)*x)-cos(N*x)-cos(x)) ...
    );
theoryBp                    = f_theoryBp(duVec,r,N);
theoryBp(isnan(theoryBp))   = (r^2)/((r-1)^2);
theoryBp_dbAbs              = db(abs(theoryBp));
theoryBp_dbAbs_norm         = theoryBp_dbAbs - max(theoryBp_dbAbs(:));
end