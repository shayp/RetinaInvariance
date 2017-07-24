% learnAndPredictAllNeuronsWithHistory
clear all;
addpath('preProcessTools');
load('globalParams');
load ('repeatStimulusTimes');    
load ('RepStimulusExtended');
load ('RepSpTimes');    
load('sortedQualityInd');
numOfNeurons = 100;
scaledRepStimulus = ShrinkRepeatStimilus(RepStimulusExtended, repeatStimulusTimes, spikesWantedSampFactor);
% for i = 1:numOfNeurons
%     neuronsInNetwork = [sortedQualityInd(i)]
%     learnAndPrediectForNetworkConfiguration(neuronsInNetwork, scaledRepStimulus)
% end
neuronsInNetwork = [33]
learnAndPrediectForNetworkConfiguration(neuronsInNetwork, scaledRepStimulus)

