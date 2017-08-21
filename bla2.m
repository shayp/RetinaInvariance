addpath('preProcessTools');
lastPeak = 0.03;
dt = 0.001;
absoluterRefractory = 0.006;
hpeaks = [absoluterRefractory + 4 * dt lastPeak];
b = -7*dt;
numOfBaseVectors = 10;
% build post spike BaseVectors
[~, ~, originalBaseVectors] = buildBaseVectorsForPostSpikeAndCoupling(numOfBaseVectors,dt,hpeaks, b, absoluterRefractory);
size(originalBaseVectors)

W = [-10 -3 -2 0.1 0 0 0 0 0 0];
% Change the resolution of the base vectors to be 80ms~
figure();
plot(originalBaseVectors);drawnow;
% load('spikes');
%  refreactoryPeriodArr =  getRefractoryPeriodForNeurons(spikes)
