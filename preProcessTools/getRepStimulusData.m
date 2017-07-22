function repSession = getRepStimulusData(networkNeurons, repeatStimulusTimes, RepSpTimes, wantedSampleFactor, RepStimulusExtended)
numOfNeurons = length(networkNeurons);
    for i = 1:numOfNeurons
        repSession(i).scaledRepSpikes  = shrinkRepeatSpikes(repeatStimulusTimes, RepSpTimes(coupledNeurons(i)).sp, wantedSampleFactor);
    end

end