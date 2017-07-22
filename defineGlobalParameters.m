function defineGlobalParameters()
stimulusFilterParamsSize = 40;
spikesWantedSampFactor = 10;
binsInSecond = 1000 / spikesWantedSampFactor;
stimulusWantedSampleFactor = spikesWantedSampFactor * 10;
numOfBaseVectors = 10;
baseVectorLength = 40;
trainFrac = 0.9;
% Define number of base vectors for post spike filter
lastPeak = 0.05;
dt = 0.001;
hpeaks = [0.001 lastPeak];
b = 0.005;
numOfRepeats = 200;
save('globalParams', 'stimulusFilterParamsSize', 'spikesWantedSampFactor', 'stimulusWantedSampleFactor', 'numOfBaseVectors',...
    'baseVectorLength', 'lastPeak','dt', 'hpeaks', 'b', 'trainFrac','binsInSecond', 'numOfRepeats');
end