function [interpSpikes, interpStimulus] = changeSpikesAndStimulusRsolution(spikes, stimulus,stimtimes, wantedSampFactor)
    endStimulus = find(stimtimes > spikes(length(spikes)));
    lastStimulus = stimtimes(endStimulus(1));
    spikesSampleVector = ismember(1:lastStimulus, spikes);
    stimulusSampleVector = zeros(lastStimulus, 1);
    for i = 1:endStimulus(1)  - 1
        stimulusSampleVector(stimtimes(i):stimtimes(i + 1)) = stimulus(i);
    end
    strech = ceil(lastStimulus / wantedSampFactor);
    interpStimulus = imresize(stimulusSampleVector, [strech 1], 'nearest');
    interpSpikes = zeros(strech,1);
    for i = 1:strech  -1
        interpSpikes(i) = sum(spikesSampleVector((i - 1) * wantedSampFactor + 1: i * wantedSampFactor));
    end
   
end