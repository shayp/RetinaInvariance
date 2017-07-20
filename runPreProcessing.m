datdir = './';  
defineGlobalParameters();
load([datdir, 'Stim']);    
load([datdir,'stimtimes']); 
load([datdir, 'SpTimes']); 
load([datdir, 'globalParams']); 

[spikes, stimulus, stimulusDesignMatrix] = BuildGeneralDataForLearning(Stim, stimtimes, SpTimes,stimulusFilterParamsSize,...
     spikesWantedSampFactor, stimulusWantedSampleFactor, baseVectorLength);

% build post spike BaseVectors
[~, ~, originalBaseVectors] = buildBaseVectorsForPostSpikeAndCoupling(numOfBaseVectors,dt,hpeaks, b);

postSpikeBaseVectors = imresize(originalBaseVectors, [baseVectorLength numOfBaseVectors]);
save('postSpikeBaseVectors', 'postSpikeBaseVectors');
