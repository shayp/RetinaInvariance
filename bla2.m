
lastPeak = 0.05;
dt = 0.001;
absoluterRefractory = 0.004;
hpeaks = [absoluterRefractory + 2 * dt lastPeak];
b = 0.001;
numOfBaseVectors = 10;
% build post spike BaseVectors
[~, ~, originalBaseVectors] = buildBaseVectorsForPostSpikeAndCoupling(numOfBaseVectors,dt,hpeaks, b, absoluterRefractory);
size(originalBaseVectors)

% Change the resolution of the base vectors to be 80ms~
figure();
plot(originalBaseVectors);drawnow;
