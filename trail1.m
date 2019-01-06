clear all;
close all;
clc;

phaseVec = linspace(0,2*pi,1000);

% 'color',rand(1,3)
figure;
hold on;
for n = 2 : 10
    plot([sin(n*phaseVec(:))./sin(phaseVec(:)) (n/2)*ones(length(phaseVec),1)], 'color',rand(1,3));
    xline(1/n);
end

legend();