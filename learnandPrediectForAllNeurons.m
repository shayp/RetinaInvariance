function learnandPrediectForAllNeurons()
    load('sortedQualityInd');
    load( 'SpTimes'); 
    load ('repeatStimulusTimes');    
    load ('RepStimulusExtended');
    load ('RepSpTimes'); 
    load('globalParams');
    scaledRepStimulus = ShrinkRepeatStimilus(RepStimulusExtended, repeatStimulusTimes, spikesWantedSampFactor);

    numOfNeurons = length(SpTimes);

    for i = 1:numOfNeurons
        networkParameters = learnModelsParameters(i);
        repSession = getRepStimulusData(networkNeurons, repeatStimulusTimes, RepSpTimes, wantedSampleFactor, RepStimulusExtended)
    end
end