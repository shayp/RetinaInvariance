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
choosedNeuron = randi(100, 1, 4);
sortedQualityInd(choosedNeuron)
couplenNeurons = [sortedQualityInd(choosedNeuron)];
%couplenNeurons = [33];
numOfNeurons = length(couplenNeurons);
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
for i = 1:numOfNeurons
    scaledRepSpikes(i).vector  = shrinkRepeatSpikes(repeatStimulusTimes, RepSpTimes(couplenNeurons(i)).sp, wantedSampleFactor);
end
%scaledRepSpikes1 = shrinkRepeatSpikes(repeatStimulusTimes, RepSpTimes(couplenNeurons(1)).sp, wantedSampleFactor);
%scaledRepSpikes2 = shrinkRepeatSpikes(repeatStimulusTimes, RepSpTimes(couplenNeurons(2)).sp, wantedSampleFactor);
scaledRepStimulus = ShrinkRepeatStimilus(RepStimulusExtended, repeatStimulusTimes, wantedSampleFactor);
%stimulusDesignMatrix = buildStimulusDesignMatrix(length(learnedSTA), scaledRepStimulus);
% repSTA = zeros(length(learnedSTA), 1);
% for i = 1:size(scaledRepSpikes1,1)
% repSTA = repSTA +  calculateSTA(stimulusDesignMatrix, scaledRepSpikes1(i,:)');
% end
% repSTA = repSTA / size(scaledRepSpikes1,1);
% figure();
% plot(1:length(learnedSTA), learnedSTA, 1:length(learnedSTA), repSTA);
%%
load('Filters.mat');
for i = 1:numOfNeurons
    cleanArray = zeros(numOfRepeats, length(scaledRepStimulus));
    neuronSim(i).simulation = cleanArray;
end
for j = 1:numOfRepeats
    response = RunGLMSimulation(numOfNeurons, scaledRepStimulus, Filters, stimulusFilterLength, couplingFilterLength, deltaT);
    %response = runSimulationPoisson(numOfNeurons, scaledRepStimulus, Filters, stimulusFilterLength, couplingFilterLength,deltaT);
    for i = 1:numOfNeurons
        neuronSim(i).simulation(j,:) =  response(i,:);
    end
end
figure();
for i = 1:numOfNeurons
[spikeRate, correlation] = CalculateCorrelatedSpikeRate(numOfRepeats, scaledRepSpikes(i).vector, neuronSim(i).simulation, 8);
    lengthOfRepeat = size(spikeRate,2);
    subplot(numOfNeurons,1,i);
    plot(2* (1:lengthOfRepeat), spikeRate(1,:), 2 * (1:lengthOfRepeat), spikeRate(2,:));
%     randPoint = randi(lengthOfRepeat - 200);
%     xlim([randPoint randPoint + 200]);
    xlim([100 500]);
    title(['R ' num2str(correlation) ' Neuron ' num2str(couplenNeurons)]);
    xlabel('Time (ms)');
    ylabel('Firing rate (spikes/sec) ');
end
