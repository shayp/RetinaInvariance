function defineGlobalParameters()
expSampling = 10000;
stimulusSpikeRatio = 2;
spikesWantedSampFactor = 10;
spikeBinSizeInms = spikesWantedSampFactor / 10;
stimulusFilterParamsSize = ceil(400 / spikeBinSizeInms / stimulusSpikeRatio);
baseVectorLength = ceil(40 / spikeBinSizeInms);
stimulusWantedSampleFactor = spikesWantedSampFactor * stimulusSpikeRatio;
binsInSecond = expSampling / spikesWantedSampFactor;
deltaT = 1 / binsInSecond;
numOfBaseVectors = 6;
stimulusFilterSizeForSimulation = stimulusFilterParamsSize * stimulusSpikeRatio;
windowSizeForFiringRate = 20 / spikeBinSizeInms;
trainFrac = 0.8;
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