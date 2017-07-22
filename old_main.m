% Load the data for GLM
%%
clear all;
datdir = './';  
load([datdir, 'Stim']);    
load([datdir,'stimtimes']); 
load([datdir, 'SpTimes']);    

%coupledNeurons = [33];
numOfNeurons = length(coupledNeurons);
wantedSampleFactor = 20;

%% Repeat stimulus
load('globalParams.mat');
load('NeuronParameters.mat');

for i = 1:numOfNeurons
    GLM_Full_NeuronParameters(i).simulation = zeros(numOfRepeats, length(scaledRepStimulus));
    GLM_Partial_NeuronParameters(i).simulation = zeros(numOfRepeats, length(scaledRepStimulus));
    LN_NeuronParameters(i).simulation = zeros(numOfRepeats, length(scaledRepStimulus));
end

for j = 1:numOfRepeats
    % Full glm simulation
    response_GLM_Full = RunSimulation(numOfNeurons, scaledRepStimulus, GLM_Full_NeuronParameters, stimulusFilterLength, couplingFilterLength, deltaT, 1);
    responseGLM_GLM_Partial = RunSimulation(numOfNeurons, scaledRepStimulus, GLM_Partial_NeuronParameters, stimulusFilterLength, couplingFilterLength, deltaT , 1);
    response_LN = RunSimulation(numOfNeurons, scaledRepStimulus, LN_NeuronParameters, stimulusFilterLength, couplingFilterLength, deltaT, 0);

    for i = 1:numOfNeurons
        NeuronParameters(i).GLMFullSimulation(j,:) =  response_GLM_Full(i,:);
        NeuronParameters(i).GLMPartialSimulation(j,:) =  responseGLM_GLM_Partial(i,:);
        NeuronParameters(i).LNSimulation(j,:) =  response_LN(i,:);
    end
end

numbrOfBins = 10;
for i = 1:numOfNeurons
    [glmFullSpikeRate, glmFullCorrelation] = CalculateCorrelatedSpikeRate(numOfRepeats, NeuronParameters(i).scaledRepSpikes, NeuronParameters(i).GLMFullSimulation, 8);
    [glmPartialSpikeRate, glmPartialCorrelation] = CalculateCorrelatedSpikeRate(numOfRepeats, NeuronParameters(i).scaledRepSpikes, NeuronParameters(i).GLMPartialSimulation, 8);
    [lnSpikeRate,lnCorrelation] = CalculateCorrelatedSpikeRate(numOfRepeats, NeuronParameters(i).scaledRepSpikes, NeuronParameters(i).LNSimulation, 8);

    NeuronParameters(i).realSpikeRate = glmFullSpikeRate(1,:);
    NeuronParameters(i).glmFullSimulatedSpikeRate = glmFullSpikeRate(2,:);
    NeuronParameters(i).glmPartialSimulatedSpikeRate = glmPartialSpikeRate(2,:);
    NeuronParameters(i).lnSimulatedSpikeRate = lnSpikeRate(2,:);

    correaltionVector = zeros(4, numbrOfBins);
    maxRealSpike = max(NeuronParameters(i).realSpikeRate) + 0.0001;
    firingRateSpace = linspace(0, maxRealSpike, numbrOfBins + 1);
    for bin = 1:numbrOfBins
        wantedIndexes = find(NeuronParameters(i).realSpikeRate >= firingRateSpace(bin) & NeuronParameters(i).realSpikeRate < firingRateSpace(bin + 1));
        sizeOfBin = length(wantedIndexes);
        correaltionVector(1, bin) = sum(NeuronParameters(i).realSpikeRate(wantedIndexes)) / sizeOfBin;
        correaltionVector(2, bin) = sum(NeuronParameters(i).glmPartialSimulatedSpikeRate(wantedIndexes)) / sizeOfBin;
        correaltionVector(3, bin) = sum(NeuronParameters(i).glmFullSimulatedSpikeRate(wantedIndexes)) / sizeOfBin;
        correaltionVector(4, bin) = sum(NeuronParameters(i).lnSimulatedSpikeRate(wantedIndexes)) / sizeOfBin;

    end
    NeuronParameters(i).correaltionVector = correaltionVector;
    NeuronParameters(i).spikeRateCorrelation = [glmPartialCorrelation glmFullCorrelation lnCorrelation]
end
save('FinalNeuronParameters.mat', 'NeuronParameters', 'GLM_Full_NeuronParameters', 'GLM_Partial_NeuronParameters', 'LN_NeuronParameters');
plotResults();
