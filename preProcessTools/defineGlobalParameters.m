function defineGlobalParameters()
expSampling = 10000;
stimulusSpikeRatio = 2;
spikesWantedSampFactor = 10;
spikeBinSizeInms = spikesWantedSampFactor / 10;
stimulusFilterParamsSize = ceil(400 / spikeBinSizeInms / stimulusSpikeRatio);
stimulusWantedSampleFactor = spikesWantedSampFactor * stimulusSpikeRatio;
binsInSecond = expSampling / spikesWantedSampFactor;
deltaT = 1 / binsInSecond;
stimulusFilterSizeForSimulation = stimulusFilterParamsSize * stimulusSpikeRatio;
windowSizeForFiringRate = 20 / spikeBinSizeInms;
trainFrac = 0.7;
% Define number of base vectors for post spike filter
lastPeak = 0.05;
numOfBaseVectors = 10;
dt = 0.001;
absoluterRefractory = 0.004;
hpeaks = [absoluterRefractory + 2 * dt lastPeak];
b = 0.001;
numOfRepeats = 10;
save('globalParams', 'stimulusFilterParamsSize', 'spikesWantedSampFactor', 'stimulusWantedSampleFactor',...
    'numOfBaseVectors', 'lastPeak','dt', 'hpeaks', 'b', 'trainFrac','binsInSecond',...
    'numOfRepeats', 'stimulusFilterSizeForSimulation', 'deltaT', 'windowSizeForFiringRate', 'stimulusSpikeRatio', 'absoluterRefractory');
end