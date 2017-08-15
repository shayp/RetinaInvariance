
lastPeak = 0.075;
dt = 0.001;
absoluterRefractory = 0.003;
ProbRefractory = absoluterRefractory;
hpeaks = [ProbRefractory + dt lastPeak];
b = 0.00001;
numOfBaseVectors = 10;
% build post spike BaseVectors
[~, ~, originalBaseVectors] = buildBaseVectorsForPostSpikeAndCoupling(numOfBaseVectors,dt,hpeaks, b, absoluterRefractory, ProbRefractory);
size(originalBaseVectors)

% Change the resolution of the base vectors to be 80ms~
figure();
plot(originalBaseVectors);drawnow;
