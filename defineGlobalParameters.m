function defineGlobalParameters()
stimulusFilterParamsSize = 20;
spikesWantedSampFactor = 20;
stimulusWantedSampleFactor = spikesWantedSampFactor * 10;
numOfBaseVectors = 10;
baseVectorLength = 40;
% Define number of base vectors for post spike filter
lastPeak = 0.05;
dt = 0.001;
hpeaks = [0.001 lastPeak];
b = 0.005;
save('globalParams', 'stimulusFilterParamsSize', 'spikesWantedSampFactor', 'stimulusWantedSampleFactor', 'numOfBaseVectors',...
    'baseVectorLength', 'lastPeak','dt', 'hpeaks', 'b');
end