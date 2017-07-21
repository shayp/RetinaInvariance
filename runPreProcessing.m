datdir = './';  
defineGlobalParameters();
load([datdir, 'Stim']);    
load([datdir,'stimtimes']); 
load([datdir, 'SpTimes']); 
load([datdir, 'globalParams']); 

% Build general data for training(Spike train, stimulus, stimulus design
% matrix)
[spikes, stimulus, stimulusDesignMatrix, postSpikeBaseVectors, spikeHistoryData] = BuildGeneralDataForLearning(Stim, stimtimes, SpTimes,stimulusFilterParamsSize,...
     spikesWantedSampFactor, stimulusWantedSampleFactor, baseVectorLength);
 
