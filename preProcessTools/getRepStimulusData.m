function repSession = getRepStimulusData(networkNeurons, repeatStimulusTimes, RepSpTimes, wantedSampleFactor)
numOfNeurons = length(networkNeurons);
    for i = 1:numOfNeurons
        repSession(i).scaledRepSpikes  = shrinkRepeatSpikes(repeatStimulusTimes, RepSpTimes(networkNeurons(i)).sp, wantedSampleFactor);
    end
end
