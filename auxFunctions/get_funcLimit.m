function [symOut,report] = get_funcLimit(symIn,limParam,limParamVal)
[num,den]   = numden(symIn);
foundExpr   = 0;
iterId      = 0;
while ~foundExpr
    numVal      = subs(num, limParam, limParamVal);
    denVal      = subs(den, limParam, limParamVal);
    if ~(denVal==0)
        foundExpr   = 1;
    else
        num     = diff(num, limParam);
        den     = diff(den, limParam);
        iterId  = iterId + 1;
    end
end
symOut     = simplify(numVal/denVal);

report.nIter = iterId;
end