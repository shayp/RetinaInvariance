% learnAndPredictAllNeuronsWithHistory
clear all;
addpath('preProcessTools');
load('globalParams');
load ('repeatStimulusTimes');    
load ('RepStimulusExtended');
load ('RepSpTimes');    
load('sortedQualityInd');
goodNeurons = [33 62 80 144 146 151 162 167 168 171 172 178 183 192 221];
numOfNeurons = 251;
scaledRepStimulus = ShrinkRepeatStimilus(RepStimulusExtended, repeatStimulusTimes, spikesWantedSampFactor);
for i = 1:numOfNeurons
    neuronsInNetwork = [sortedQualityInd(i)];
    learnAndPrediectForNetworkConfiguration(neuronsInNetwork, scaledRepStimulus)
end
