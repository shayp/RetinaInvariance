function [spikes, stimulus, stimulusDesignMatrix, postSpikeBaseVectors] = BuildGeneralDataForLearning(Stim, stimtimes, SpTimes,stimulusFilterParamsSize,...
    spikesWantedSampFactor, stimulusWantedSampleFactor, numOfBaseVectors, baseVectorLength)


stimulus = changeStimulusRsolution(Stim,stimtimes, stimulusWantedSampleFactor);
save('stimulus', 'stimulus');

lastIndex = length(stimulus) * stimulusWantedSampleFactor;

stimulusDesignMatrix = buildStimulusDesignMatrix(stimulusFilterParamsSize, stimulus);
save('stimulusDesignMatrix', 'stimulusDesignMatrix');

spikes = changeSpikeRsolution(SpTimes, lastIndex, spikesWantedSampFactor);
save('spikes', 'spikes');

end