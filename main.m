% Load the data for GLM
%%
clear all;
datdir = './';  
load([datdir, 'Stim']);    
load([datdir,'stimtimes']); 
load([datdir, 'SpTimes']); 
load ('repeatStimulusTimes');    
load ('RepStimulusExtended');
load ('RepSpTimes');    
load('sortedQualityInd');
ncells = length(SpTimes);
numOfRepeats = 200;
choosedNeuron = randi(6);
sortedQualityInd(choosedNeuron)
couplenNeurons = [sortedQualityInd(choosedNeuron)];
%couplenNeurons = [33];
numOfNeurons = 1;
wantedSampleFactor = 20;
%%
for i = 1:length(couplenNeurons)
[scaledStimulus, couplingFilters, learnedSTA, deltaT, meanFiringRate] = runGLM(couplenNeurons(i), Stim, stimtimes, SpTimes, couplenNeurons);
stimulusFilterLength = length(learnedSTA);
couplingFilterLength = size(couplingFilters,2);
Filters(i).StimulusFilter = learnedSTA;
Filters(i).couplingFilters = couplingFilters;
Filters(i).meanFiringRate = meanFiringRate;
end
save('Filters.mat', 'Filters', 'scaledStimulus', 'stimulusFilterLength', 'couplingFilterLength', 'deltaT');
%% Repeat stimulus
scaledRepSpikes1 = shrinkRepeatSpikes(repeatStimulusTimes, RepSpTimes(couplenNeurons(1)).sp, wantedSampleFactor);
%scaledRepSpikes2 = shrinkRepeatSpikes(repeatStimulusTimes, RepSpTimes(couplenNeurons(2)).sp, wantedSampleFactor);
scaledRepStimulus = ShrinkRepeatStimilus(RepStimulusExtended, repeatStimulusTimes, wantedSampleFactor);
stimulusDesignMatrix = buildStimulusDesignMatrix(length(learnedSTA), scaledRepStimulus);
repSTA = zeros(length(learnedSTA), 1);
for i = 1:size(scaledRepSpikes1,1)
repSTA = repSTA +  calculateSTA(stimulusDesignMatrix, scaledRepSpikes1(i,:)');
end
repSTA = repSTA / size(scaledRepSpikes1,1);
figure();
plot(1:length(learnedSTA), learnedSTA, 1:length(learnedSTA), repSTA);
%%
neuron1Sim = zeros(numOfRepeats, length(scaledRepStimulus));
neuron2Sim = zeros(numOfRepeats, length(scaledRepStimulus));
load('Filters.mat');

for j = 1:numOfRepeats
    response = RunGLMSimulation(numOfNeurons, scaledRepStimulus, Filters, stimulusFilterLength, couplingFilterLength, deltaT);
    %response = runSimulationPoisson(numOfNeurons, scaledRepStimulus, Filters, stimulusFilterLength, couplingFilterLength,deltaT);
    neuron1Sim(j,:) =  response(1,:);
    %neuron2Sim(j,:) =  response(2,:);
end

[spikeRate, correlation] = CalculateCorrelatedSpikeRate(numOfRepeats, scaledRepSpikes1, neuron1Sim, 8);
lengthOfRepeat = size(spikeRate,2);
figure();
plot(2* (1:lengthOfRepeat), spikeRate(1,:), 2 * (1:lengthOfRepeat), spikeRate(2,:));
randPoint = randi(lengthOfRepeat - 800);
xlim([randPoint randPoint + 800]);
%xlim([1 500]);
title(['R ' num2str(correlation)]);
xlabel('Time (ms)');
ylabel('Firing rate ');
