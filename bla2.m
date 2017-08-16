
lastPeak = 0.025;
dt = 0.001;
absoluterRefractory = 0.003;
hpeaks = [absoluterRefractory + 2 * dt lastPeak];
b = 0.00001;
numOfBaseVectors = 10;
% build post spike BaseVectors
[~, ~, originalBaseVectors] = buildBaseVectorsForPostSpikeAndCoupling(numOfBaseVectors,dt,hpeaks, b, absoluterRefractory);
size(originalBaseVectors)

% Change the resolution of the base vectors to be 80ms~
figure();
plot(originalBaseVectors);drawnow;

% load('spikes');
%  refreactoryPeriodArr =  getRefractoryPeriodForNeurons(spikes)
