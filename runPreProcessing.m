load( 'NonRepeatStim');    
load('NonRepeatstimtimes'); 
load('SpTimes'); 
addpath('preProcessTools');
defineGlobalParameters();
load('globalParams'); 

% Build general data for training(Spike train, stimulus, stimulus design
% matrix)
[spikes, stimulus, stimulusDesignMatrix, postSpikeBaseVectors, spikeHistoryData] = BuildGeneralDataForLearning(Stim, stimtimes, SpTimes,stimulusFilterParamsSize,...
     spikesWantedSampFactor, stimulusWantedSampleFactor, baseVectorLength);
 
