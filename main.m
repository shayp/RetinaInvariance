% Load the data for GLM
%%
datdir = './';  
load([datdir, 'Stim']);    
load([datdir,'stimtimes']); 
load([datdir, 'SpTimes']); 
load ('repeatStimulusTimes');    
load ('RepStimulusExtended');
load ('RepSpTimes');    
ncells = length(SpTimes);
numOfRepeats = length(repeatStimulusTimes);
couplenNeurons = [1 2];
numOfNeurons = 2;
wantedSampleFactor = 20;
%%
[scaledStimulus, couplingFilters, learnedSTA, meanFiringRate] = runGLM(couplenNeurons(1), Stim, stimtimes, SpTimes, couplenNeurons);
stimulusFilterLength = length(learnedSTA);
couplingFilterLength = size(couplingFilters,2);
Filters(1).StimulusFilter = learnedSTA;
Filters(1).couplingFilters = couplingFilters;
Filters(1).meanFiringRate = meanFiringRate;


[scaledStimulus, couplingFilters, learnedSTA, deltaT, meanFiringRate] = runGLM(couplenNeurons(2), Stim, stimtimes, SpTimes, couplenNeurons);
Filters(2).StimulusFilter = learnedSTA;
Filters(2).couplingFilters = couplingFilters;
Filters(2).meanFiringRate = meanFiringRate;

save('Filters.mat', 'Filters', 'scaledStimulus', 'stimulusFilterLength', 'couplingFilterLength', 'deltaT');
%% Repeat stimulus
scaledRepSpikes1 = shrinkRepeatSpikes(repeatStimulusTimes, RepSpTimes(couplenNeurons(1)).sp, wantedSampleFactor);
scaledRepSpikes2 = shrinkRepeatSpikes(repeatStimulusTimes, RepSpTimes(couplenNeurons(2)).sp, wantedSampleFactor);
scaledRepStimulus = ShrinkRepeatStimilus(RepStimulusExtended, repeatStimulusTimes, wantedSampleFactor);


%%

neuron1Sim = zeros(numOfRepeats, length(scaledRepStimulus));
neuron2Sim = zeros(numOfRepeats, length(scaledRepStimulus));
load('Filters.mat');
for j = 1:numOfRepeats
    response = RunGLMSimulation(numOfNeurons, scaledRepStimulus, Filters, stimulusFilterLength, couplingFilterLength, deltaT);
    neuron1Sim(j,:) =  response(1,:);
    neuron2Sim(j,:) =  response(2,:);
end
spikeRate = CalculateCorrelatedSpikeRate(numOfRepeats, scaledRepSpikes2, neuron2Sim, 15);
lengthOfRepeat = size(spikeRate,2);
figure();
plot(1:lengthOfRepeat, spikeRate(1,:), 1:lengthOfRepeat, spikeRate(2,:));
xlim([400 600]);