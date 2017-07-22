function [stimulus] = changeStimulusResolution(stimulusVector,stimtimes, wantedSampFactor)
    lastStimulus = stimtimes(end);
    stimulusSampleVector = zeros(lastStimulus, 1);
    for i = 1:length(stimulusVector) - 1
        stimulusSampleVector(stimtimes(i):stimtimes(i + 1) - 1) = stimulusVector(i);
    end
    strech = ceil(lastStimulus / wantedSampFactor);
    stimulus = imresize(stimulusSampleVector, [strech 1], 'nearest');
end