function scaledRepStimulus = ShrinkRepeatStimilus(stimulus, stimTimes, wantedSampleFactor)
    currentLength = stimTimes(2) - stimTimes(1);
    LongStim = stimulus(1,stimTimes(1) + 1: stimTimes(2));
    wantedLength = ceil(currentLength / wantedSampleFactor)
    scaledRepStimulus = zeros(wantedLength, 1);

    for i = 1:wantedLength - 1
        scaledRepStimulus(i) = mean(LongStim((i - 1) * wantedSampleFactor + 1: i * wantedSampleFactor));
    end
        scaledRepStimulus(i) = mean(LongStim((wantedLength - 1) * wantedSampleFactor + 1 :end));
end