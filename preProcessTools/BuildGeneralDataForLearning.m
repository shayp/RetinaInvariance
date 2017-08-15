function [spikes, stimulus, stimulusDesignMatrix, postSpikeBaseVectors, spikeHistoryData] = BuildGeneralDataForLearning(Stim, stimtimes, SpTimes,stimulusFilterParamsSize,...
    spikesWantedSampFactor, stimulusWantedSampleFactor, numOfBaseVectors)
load('globalParams'); 

stimulus = changeStimulusResolution(Stim(1:end - 1),stimtimes(1:end -1), stimulusWantedSampleFactor);




lastIndex = length(stimulus) * stimulusWantedSampleFactor;

stimulusDesignMatrix = buildStimulusDesignMatrix(stimulusFilterParamsSize, stimulus);
save('stimulusDesignMatrix', 'stimulusDesignMatrix');

spikes = changeSpikeResolution(SpTimes, lastIndex, spikesWantedSampFactor);
spikesLength = length(spikes(1).data);
save('spikes', 'spikes');

% Get the fine stimulus resolution(Just for calculating the STA)
fineStimulus = changeStimulusResolution(Stim,stimtimes, spikesWantedSampFactor);

fineStimulus = fineStimulus(1:spikesLength);
save('stimulus', 'stimulus', 'fineStimulus');

% Calculate the fine stimulus DesignMatrix(For STA only)
fineStimulusDesignMatrix = buildStimulusDesignMatrix(stimulusFilterParamsSize * (stimulusWantedSampleFactor / spikesWantedSampFactor), fineStimulus);

[stimlusFilters, fineStimulusFilters] = getStimuliusFilterForAllCells(spikes, fineStimulusDesignMatrix, stimulusFilterParamsSize);
save('stimlusFilters', 'stimlusFilters', 'fineStimulusFilters');

% build post spike BaseVectors
[~, ~, postSpikeBaseVectors] = buildBaseVectorsForPostSpikeAndCoupling(numOfBaseVectors,dt,hpeaks, b, absoluterRefractory);

figure();
plot(postSpikeBaseVectors);drawnow;
save('postSpikeBaseVectors', 'postSpikeBaseVectors','numOfBaseVectors');

spikeHistoryData = getSpikeHistoryDataForNeurons(spikes, numOfBaseVectors, postSpikeBaseVectors);
save('spikeHistoryData', 'spikeHistoryData', '-v7.3');

end