function [sigData] = generateSignalData(cfgIn)
%% configure
defautCfg.duration_SAMPLES      = 1024;
defautCfg.baseFreq_relative     = 0.01;
defautCfg.bandwidth_relative    = 0.001;

cfgFields = fieldnames(defautCfg);

for cfgFieldID = 1 : numel(cfgFields)
    curCfgField     = cfgFields{cfgFieldID};
    try
        cmdString   = [curCfgField '=cfgIn.(''' curCfgField ''');'];
        eval(cmdString);
    catch
        cmdString   = [curCfgField '=defautCfg.(''' curCfgField ''');'];
        eval(cmdString);
    end
end

%% generate data
relativeFreqVec     = linspace(...
    baseFreq_relative-0.5*bandwidth_relative, ...
    baseFreq_relative+0.5*bandwidth_relative, ...
    duration_SAMPLES);

phaseIncVec     = 2*pi*reshape(relativeFreqVec,1,[]);
phaseVec        = cumsum([0 phaseIncVec(1:end-1)]);
sigData         = exp(1i*phaseVec);

if false
    %% DEBUG
    sigFFT      = db(abs(fftshift(fft(sigData))));
    figure;
    subplot(2,1,1);
    plot([real(sigData(:)) imag(sigData(:))]);
    subplot(2,1,2);
    freqVec     = linspace(-0.5,0.5,duration_SAMPLES);
    plot(freqVec, sigFFT);
end
end