function [] = fisherInformationMatrixCalc()
%% configure
cVal            = 3e8;
fVal            = 10e9;
RVal            = 10e3;
noisePowerDB    = -60;
psVal           = pi/2;
lambdaVal       = cVal/fVal;
if true
    %% N
    init_N      = 2;
    final_N     = 200;
    hop_N       = 5;
    %% r
    nVal_r      = 30;
    init_r      = 0.6;
    final_r     = 0.99;
    %% DU
    nVal_DU     = 40;
    init_DU     = 0;
    final_DU    = 1;
    %% DT
    nVal_DT     = 30;
    init_DT     = 0;
    final_DT    = .04*lambdaVal/cVal;
    %% w
    nVal_w     = 30;
    init_w     = 2*pi*fVal;
    final_w    = 2*pi*1.1*fVal;
end
N = init_N;

%% symbolics
syms D DU ts DT c w r ps p s positive;
%{
D   : inter element spacing
DU  : electrical phase diff
ts  : expected propagation delay
DT  : propagation delay diff
c   : propagation velocity
w   : signal radial frequency
r   : array misalignemnt factor
ps  : expected DOA
s   : noise variance
%}

%% basic terms
tauP            = D*cos(ps)/c;
us              = 2*pi*w*tauP;
u               = us + DU;
sensorIdVec     = 0:(N-1);
ds              = reshape(exp(-1i*sensorIdVec*us),[],1);
dsT             = transpose(ds);
d               = reshape(exp(-1i*sensorIdVec*u),[],1);
dT              = transpose(d);
t               = ts + DT;

%% dual convetional beamformer (DCBF)
DCBF_a      = (r/N)*conj(ds)*exp(-1i*ts);
DCBF_aT     = transpose(DCBF_a);
DCBF_aH     = conj(DCBF_aT);
DCBF_b      = (4*c^2*ts^2*r/N)*conj(ds)*exp(-1i*ts);
DCBF_bT     = transpose(DCBF_b);
DCBF_bH     = conj(DCBF_bT);

%% FIM (singal freq)
%{
pp  - theta,    theta
pt  - theta,    tau
tp  - tau,      theta
tt  - tua,      tau
%}
a   = sym('a',[N,1]);
aT  = transpose(a);
aH  = conj(aT);
b   = sym('b',[N,1]);
bT  = transpose(b);
bH  = conj(bT);
S   = (1/(2*pi*s));
if true
    %% basic terms
    A       = -1i*w*tauP*diag(sensorIdVec);
    B       = d*dT*A - A*d*dT;
    jDen    = (4*c^2^t^2-bT*d*exp(-1i*w*t))^2;
    %% pp
    ppTerm1     = aT*((4*c^2*t^2)*A*d + B*b*exp(-1i*w*t));
    ppTerm2     = ppTerm1/jDen;
    pp          = S*ppTerm2*conj(ppTerm2);
    %% tt
    ttTerm1     = w^2*aT*d*conj(aT*d);
    ttTerm2     = ttTerm1/jDen;
    tt          = S*ttTerm2*conj(ttTerm2);
    %% pt , tp
    ptTerm1     = real(1i*w*aT*(A*d+B*b*exp(-1i*w*t))*aH*conj(d));
    pt          = S*ptTerm1/jDen;
    tp          = pt;
    %% J
    J   = [pp pt ; tp tt];
    %% CRLB
    if true
        %% calculation
        CRLB    = inv(J);
        ppCRLB  = CRLB(1,1);
        ttCRLB  = CRLB(2,2);
        %% derivatives
        ppCRLB_aDiff = sym(zeros(N,1));
        ppCRLB_bDiff = sym(zeros(N,1));
        ttCRLB_aDiff = sym(zeros(N,1));
        ttCRLB_bDiff = sym(zeros(N,1));
        for diffID = 1:N
            ppCRLB_aDiff(diffID)    = diff(ppCRLB,a(diffID));
            ppCRLB_bDiff(diffID)    = diff(ppCRLB,b(diffID));
            ttCRLB_aDiff(diffID)    = diff(ttCRLB,a(diffID));
            ttCRLB_bDiff(diffID)    = diff(ttCRLB,b(diffID));
        end
        ppCRLB_diff = [ppCRLB_aDiff ; ppCRLB_bDiff];
        ttCRLB_diff = [ttCRLB_aDiff ; ttCRLB_bDiff];
    end
end

%% simulation
if true
    %% assign the DCBF
    ppCRLB_fullPrm              = subs(ppCRLB,[a(:) ; b(:)], [DCBF_a(:) ; DCBF_b(:)]);
    ttCRLB_fullPrm              = subs(ttCRLB,[a(:) ; b(:)], [DCBF_a(:) ; DCBF_b(:)]);
    ppCRLB_diff_DCBF_fullPrm    = subs(ppCRLB_diff,[a(:) ; b(:)], [DCBF_a(:) ; DCBF_b(:)]);
    ttCRLB_diff_DCBF_fullPrm    = subs(ttCRLB_diff,[a(:) ; b(:)], [DCBF_a(:) ; DCBF_b(:)]);
    %% generate input
    rVec    = linspace(init_r,final_r,nVal_r);
    DUVec   = linspace(init_DU,final_DU,nVal_DU);
    DTVec   = linspace(init_DT,final_DT,nVal_DT);
    wVec   = linspace(init_w,final_w,nVal_w);
    %% sweep parameters
    DVal    = lambdaVal/2;
    sVal    = 10^(noisePowerDB/20);
    tsVal   = RVal/cVal;
    if true
        %% ppCRLB
        if true
            ppCRLB_DCBF    = subs(ppCRLB_fullPrm,[D, c, s, ts],[DVal, cVal, sVal, tsVal]);
            try
                lim_ppCRLB_diff_DCBF    = subs(ppCRLB,ps,psVal);
            catch
                %% steering to psVal using L'Hospital
                [num,den]   = numden(ppCRLB_DCBF);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, ps, psVal);
                        denVal      = subs(den, ps, psVal);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = diff(num, ps);
                            den         = diff(den, ps);
                        end
                    catch
                        num         = diff(num, ps);
                        den         = diff(den, ps);
                    end
                end
                curExpr                 = numVal/denVal;
                lim_ppCRLB_diff_DCBF    = curExpr
            end
            for cur_r = rVec
                curExpr1    = subs(lim_ppCRLB_diff_DCBF, r, cur_r);
                for curw = wVec
                    curExpr2     = subs(curExpr1, w, curw);
                    for curDT = DTVec
                        curExpr3     = subs(curExpr2, DT, curDT);
                        for curDU = DUVec
                            finalExpr     = eval(subs(curExpr3, DU, curDU));
                            log(abs(finalExpr))
                        end
                    end
                end
            end
        end
        %% ppCRLB_diff_DCBF
        if false
            ppCRLB_diff_DCBF    = subs(ppCRLB_diff_DCBF_fullPrm,[D, c, s, ts],[DVal, cVal, sVal, tsVal]);
            try
                lim_ppCRLB_diff_DCBF    = subs(ppCRLB_diff_DCBF,ps,psVal);
            catch
                %% steering to psVal using L'Hospital
                [num,den]   = numden(ppCRLB_diff_DCBF);
                foundExpr   = 0;
                while ~foundExpr
                    try
                        numVal      = subs(num, ps, psVal);
                        denVal      = subs(den, ps, psVal);
                        if ~(denVal==0)
                            foundExpr   = 1;
                        else
                            num         = diff(num, ps);
                            den         = diff(den, ps);
                        end
                    catch
                        num         = diff(num, ps);
                        den         = diff(den, ps);
                    end
                end
                curExpr                 = numVal/denVal;
                lim_ppCRLB_diff_DCBF    = curExpr
            end
            for cur_r = rVec
                curExpr1    = subs(lim_ppCRLB_diff_DCBF, r, cur_r);
                for curw = wVec
                    curExpr2     = subs(curExpr1, w, curw);
                    for curDT = DTVec
                        curExpr3     = subs(curExpr2, DT, curDT);
                        for curDU = DUVec
                            finalExpr     = eval(subs(curExpr3, DU, curDU));
                            log(abs(finalExpr))
                        end
                    end
                end
            end
        end
    end
end