function scaledRepStimulus = ShrinkRepeatStimilus(stimulus, stimTimes, wantedSampleFactor)
    currentLength = stimTimes(2) - stimTimes(1);
    LongStim = stimulus( stimTimes(1) + 1: stimTimes(2));
    wantedLength = ceil(currentLength / wantedSampleFactor);
    scaledRepStimulus = imresize(LongStim, [wantedLength 1], 'nearest');
end