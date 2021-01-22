function [yResample] = f_resample(x,y,xResample)
ySliced = y;
xSliced = x;
yResample = interp1(xSliced(:),ySliced(:),xResample,'spline');    
end