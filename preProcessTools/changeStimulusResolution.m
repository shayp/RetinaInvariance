function [stimulus] = changeStimulusResolution(stimulusVector,stimtimes, wantedSampFactor)
    lastStimulus = stimtimes(end);
    stimulusSampleVector = zeros(lastStimulus, 1);
    for i = 1:length(stimulusVector) - 1
        stimulusSampleVector(stimtimes(i):stimtimes(i + 1) - 1) = stimulusVector(i);
    end
    strech = ceil(lastStimulus / wantedSampFactor);
    
    stimulus = zeros(strech, 1);

    for i = 1:strech - 1
        stimulus(i) = mode(stimulusSampleVector((i - 1) * wantedSampFactor + 1: i * wantedSampFactor));
    end
    
    stimulus(strech) = mode(stimulusSampleVector((strech - 1) * wantedSampFactor + 1 :end));
end