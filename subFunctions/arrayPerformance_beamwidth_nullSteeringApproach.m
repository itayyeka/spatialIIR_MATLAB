clear all;
close all;
clc;

syms x;
syms N;
syms r;
syms t;

origFun     = ...
    ((N*sin(x)/sin(N*x))^2) ...
    -2*r*cos(t-((N-1)*x/2))*N*sin(x)/sin(N*x) ...
    ;

if true
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
    
    lim_d_dt2 = simplify((numDiffVal/denDiffVal));
end