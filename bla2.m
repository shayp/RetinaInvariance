addpath('preProcessTools');
lastPeak = 0.07;
dt = 0.001;
absoluterRefractory = 0.002;
hpeaks = [absoluterRefractory + dt lastPeak];
b = dt * 5;
numOfBaseVectors = 6;
% build post spike BaseVectors
[~, ~, originalBaseVectors] = buildBaseVectorsForPostSpikeAndCoupling(numOfBaseVectors,dt,hpeaks, b, absoluterRefractory);
size(originalBaseVectors)

W = [-10 -3 -2 0.1 0 0 0 0 0 0];
% Change the resolution of the base vectors to be 80ms~
figure();
originalBaseVectors(end - 9:end,:)
plot(originalBaseVectors);drawnow;
% load('spikes');
%  refreactoryPeriodArr =  getRefractoryPeriodForNeurons(spikes)
