function defineGlobalParameters()
expSampling = 10000;
stimulusFilterParamsSize = 40;
stimulusSpikeRatio = 10;
spikesWantedSampFactor = 10;
stimulusWantedSampleFactor = spikesWantedSampFactor * stimulusSpikeRatio;
binsInSecond = expSampling / spikesWantedSampFactor;
deltaT = 1 / binsInSecond;
numOfBaseVectors = 5;
baseVectorLength = 25;
stimulusFilterSizeForSimulation = stimulusFilterParamsSize * stimulusSpikeRatio;
windowSizeForFiringRate = 16;
trainFrac = 0.5;
% Define number of base vectors for post spike filter
lastPeak = 0.05;
dt = 0.001;
hpeaks = [0.001 lastPeak];
b = 0.005;
numOfRepeats = 200;
save('globalParams', 'stimulusFilterParamsSize', 'spikesWantedSampFactor', 'stimulusWantedSampleFactor',...
    'numOfBaseVectors','baseVectorLength', 'lastPeak','dt', 'hpeaks', 'b', 'trainFrac','binsInSecond',...
    'numOfRepeats', 'stimulusFilterSizeForSimulation', 'deltaT', 'windowSizeForFiringRate', 'stimulusSpikeRatio');
end