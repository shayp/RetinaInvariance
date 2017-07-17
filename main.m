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
choosedNeuron = randi(100, 1, 1);
sortedQualityInd(choosedNeuron)
coupledNeurons = [sortedQualityInd(choosedNeuron)];
%coupledNeurons = [33];
numOfNeurons = length(coupledNeurons);
wantedSampleFactor = 20;
%%
for i = 1:length(coupledNeurons)
[scaledStimulus, couplingFilters, learnedSTA, deltaT, meanFiringRate, realSTA] = runGLM(coupledNeurons(i), Stim, stimtimes, SpTimes, coupledNeurons);
stimulusFilterLength = length(learnedSTA);
couplingFilterLength = size(couplingFilters,2);
NeuronParameters(i).neuronIndex = coupledNeurons(i);
NeuronParameters(i).coupledNeurons = coupledNeurons;
NeuronParameters(i).expStimulusFilter = realSTA;
NeuronParameters(i).StimulusFilter = learnedSTA;
NeuronParameters(i).couplingFilters = couplingFilters;
NeuronParameters(i).meanFiringRate = meanFiringRate;
end
save('NeuronParameters.mat', 'NeuronParameters');
save('globalParams.mat','scaledStimulus', 'stimulusFilterLength', 'couplingFilterLength', 'deltaT', 'numOfNeurons');
%% Repeat stimulus
load('globalParams.mat');
load('NeuronParameters.mat');
for i = 1:numOfNeurons
    NeuronParameters(i).scaledRepSpikes  = shrinkRepeatSpikes(repeatStimulusTimes, RepSpTimes(coupledNeurons(i)).sp, wantedSampleFactor);
end

scaledRepStimulus = ShrinkRepeatStimilus(RepStimulusExtended, repeatStimulusTimes, wantedSampleFactor);

for i = 1:numOfNeurons
    NeuronParameters(i).simulation = zeros(numOfRepeats, length(scaledRepStimulus));
end

for j = 1:numOfRepeats
    response = RunGLMSimulation(numOfNeurons, scaledRepStimulus, NeuronParameters, stimulusFilterLength, couplingFilterLength, deltaT);
    for i = 1:numOfNeurons
        NeuronParameters(i).simulation(j,:) =  response(i,:);
    end
end
numbrOfBins = 10;
for i = 1:numOfNeurons
    [spikeRate, correlation] = CalculateCorrelatedSpikeRate(numOfRepeats, NeuronParameters(i).scaledRepSpikes, NeuronParameters(i).simulation, 8);
    NeuronParameters(i).realSpikeRate = spikeRate(1,:);
    NeuronParameters(i).simulatedSpikeRate = spikeRate(2,:);
    correaltionVector = zeros(2, numbrOfBins);
    maxRealSpike = max(spikeRate(1,:)) + 0.0001;
    firingRateSpace = linspace(0, maxRealSpike, numbrOfBins + 1);
    for bin = 1:numbrOfBins
        wantedIndexes = find(spikeRate(1,:) >= firingRateSpace(bin) & spikeRate(1,:) < firingRateSpace(bin + 1));
        sizeOfBin = length(wantedIndexes);
        correaltionVector(1, bin) = sum(spikeRate(1,wantedIndexes)) / sizeOfBin;
        correaltionVector(2, bin) = sum(spikeRate(2,wantedIndexes)) / sizeOfBin;
    end
    NeuronParameters(i).correaltionVector = correaltionVector;
    NeuronParameters(i).spikeRateCorrelation = correlation;
end
save('FinalNeuronParameters.mat', 'NeuronParameters');
plotResults();
