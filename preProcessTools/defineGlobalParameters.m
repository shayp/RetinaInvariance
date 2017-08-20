function defineGlobalParameters()
expSampling = 10000;
stimulusSpikeRatio = 5;
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
lastPeakHistory =0.07;
lastPeakCoupling = 0.025;
numOfBaseVectorsHistory = 10;
numOfBaseVectorsCoupling = 4;
dt = deltaT;
bForHistory = -dt;
bForCoupling = dt * 5;
numOfRepeats = 200;

save('globalParams', 'stimulusFilterParamsSize', 'spikesWantedSampFactor', 'stimulusWantedSampleFactor',...
    'numOfBaseVectorsHistory', 'numOfBaseVectorsCoupling', 'lastPeakHistory', 'lastPeakCoupling', 'dt', 'bForHistory', 'bForCoupling', 'trainFrac','binsInSecond',...
    'numOfRepeats', 'stimulusFilterSizeForSimulation', 'deltaT', 'windowSizeForFiringRate', 'stimulusSpikeRatio');
end