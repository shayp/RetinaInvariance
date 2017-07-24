function [spikes, stimulus, stimulusDesignMatrix, postSpikeBaseVectors, spikeHistoryData] = BuildGeneralDataForLearning(Stim, stimtimes, SpTimes,stimulusFilterParamsSize,...
    spikesWantedSampFactor, stimulusWantedSampleFactor, numOfBaseVectors, baseVectorLength)
load('globalParams'); 

stimulus = changeStimulusResolution(Stim(1:end - 1),stimtimes(1:end -1), stimulusWantedSampleFactor);
save('stimulus', 'stimulus');

% Get the fine stimulus resolution(Just for calculating the STA)
fineStimulus = changeStimulusResolution(Stim,stimtimes, spikesWantedSampFactor);

lastIndex = length(stimulus) * stimulusWantedSampleFactor;

stimulusDesignMatrix = buildStimulusDesignMatrix(stimulusFilterParamsSize, stimulus);
save('stimulusDesignMatrix', 'stimulusDesignMatrix');

spikes = changeSpikeResolution(SpTimes, lastIndex, spikesWantedSampFactor);
spikesLength = length(spikes(1).data);
save('spikes', 'spikes');

% Calculate the fine stimulus DesignMatrix(For STA only)
fineStimulusDesignMatrix = buildStimulusDesignMatrix(stimulusFilterParamsSize * (stimulusWantedSampleFactor / spikesWantedSampFactor), fineStimulus(1:spikesLength));

[stimlusFilters, fineStimulusFilters] = getStimuliusFilterForAllCells(spikes, fineStimulusDesignMatrix, stimulusFilterParamsSize);
save('stimlusFilters', 'stimlusFilters', 'fineStimulusFilters');

% build post spike BaseVectors
[~, ~, originalBaseVectors] = buildBaseVectorsForPostSpikeAndCoupling(numOfBaseVectors,dt,hpeaks, b);

% Change the resolution of the base vectors to be 80ms~
postSpikeBaseVectors = imresize(originalBaseVectors, [baseVectorLength numOfBaseVectors]);
save('postSpikeBaseVectors', 'postSpikeBaseVectors','numOfBaseVectors');

spikeHistoryData = getSpikeHistoryDataForNeurons(spikes, numOfBaseVectors, postSpikeBaseVectors);
save('spikeHistoryData', 'spikeHistoryData', '-v7.3');

end