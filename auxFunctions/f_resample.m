function [yResample] = f_resample(x,y,xResample)
if length(x)>length(y)
    ySliced     = y;
    xSliced     = x(end-length(y)+1:end);
else
    ySliced     = y(end-length(x)+1:end);
    xSliced     = x;
end

yResample                   = zeros(size(xResample));
if ~isempty(xSliced)
    valIdSampleId               = logical(double(xResample>=min(xSliced)) .* double(xResample<=max(xSliced)));
    yResample(valIdSampleId)    = interp1(xSliced(:),ySliced(:),xResample(valIdSampleId),'spline');
    if false
        %% DEBUG
        figure;
        plot(xSliced,real(ySliced),'-*');
        hold on;
        plot(xResample,real(yResample));
        N    = size(xResample,2);
        legendStr   = cellfun(@(sensorId) ['sensor #' num2str(sensorId)], num2cell(0:(N-1)), 'UniformOutput', false);
        legend(['original' ; legendStr(:)]);
        try
            xlim([xResample(1) xResample(100)]);
        catch
        end
        close all;
    end
end
end