function learnAndPrediectForNetworkConfiguration(neuronsInNetwork, scaledRepStimulus)
    load('sortedQualityInd');
    load( 'SpTimes'); 
    load ('repeatStimulusTimes');    
    load ('RepStimulusExtended');
    load ('RepSpTimes'); 
    load('globalParams');
    

    numOfNeurons = length(neuronsInNetwork);
    learnedParameters = learnModelsParameters(neuronsInNetwork);
    NeuronParameters = getRepStimulusData(neuronsInNetwork, repeatStimulusTimes, RepSpTimes, wantedSampleFactor, RepStimulusExtended);
    
    for i = 1:numOfNeurons
        glmFullParams(i).neuronIndex = learnedParameters(i).neuronNumber;
        glmPartialParams(i).neuronIndex = learnedParameters(i).neuronNumber;
        lnParams(i).neuronIndex = learnedParameters(i).neuronNumber;

        glmFullParams(i).stimulusFilter = learnedParameters(i).fullGLMParams.StimulusFilter;
        glmPartialParams(i).stimulusFilter = learnedParameters(i).partialGLMParams.StimulusFilter;
        lnParams(i).stimulusFilter = learnedParameters(i).lnOptParams.StimulusFilter;

        glmFullParams(i).meanFiringRate = learnedParameters(i).fullGLMParams.meanFiringRate;
        glmPartialParams(i).meanFiringRate = learnedParameters(i).partialGLMParams.meanFiringRate;
        lnParams(i).meanFiringRate = learnedParameters(i).lnOptParams.meanFiringRate;

        glmFullParams(i).couplingFilters = learnedParameters(i).fullGLMParams.couplingFilters;
        glmPartialParams(i).couplingFilters = learnedParameters(i).partialGLMParams.couplingFilters;
    
        glmFullParams(i).simulation = zeros(numOfRepeats, length(scaledRepStimulus));
        glmPartialParams(i).simulation = zeros(numOfRepeats, length(scaledRepStimulus));
        lnParams(i).simulation = zeros(numOfRepeats, length(scaledRepStimulus));
    end

    for j = 1:numOfRepeats
        % Full glm simulation
        response_GLM_Full = RunSimulation(numOfNeurons, scaledRepStimulus, glmFullParams, stimulusFilterSizeForSimulation, baseVectorLength, deltaT, 1);
        responseGLM_GLM_Partial = RunSimulation(numOfNeurons, scaledRepStimulus, glmPartialParams, stimulusFilterSizeForSimulation, baseVectorLength, deltaT , 1);
        response_LN = RunSimulation(numOfNeurons, scaledRepStimulus, lnParams, stimulusFilterSizeForSimulation, baseVectorLength, deltaT, 0);
    end
    
    for i = 1:numOfNeurons
        NeuronParameters(i).GLMFullSimulation(j,:) =  response_GLM_Full(i,:);
        NeuronParameters(i).GLMPartialSimulation(j,:) =  responseGLM_GLM_Partial(i,:);
        NeuronParameters(i).LNSimulation(j,:) =  response_LN(i,:);
    end

    numbrOfBins = 10;
    for i = 1:numOfNeurons
        [glmFullSpikeRate, glmFullCorrelation] = CalculateCorrelatedSpikeRate(numOfRepeats, NeuronParameters(i).scaledRepSpikes, NeuronParameters(i).GLMFullSimulation, windowSizeForFiringRate);
        [glmPartialSpikeRate, glmPartialCorrelation] = CalculateCorrelatedSpikeRate(numOfRepeats, NeuronParameters(i).scaledRepSpikes, NeuronParameters(i).GLMPartialSimulation, windowSizeForFiringRate);
        [lnSpikeRate,lnCorrelation] = CalculateCorrelatedSpikeRate(numOfRepeats, NeuronParameters(i).scaledRepSpikes, NeuronParameters(i).LNSimulation, windowSizeForFiringRate);

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
        NeuronParameters(i).spikeRateCorrelation = [glmPartialCorrelation glmFullCorrelation lnCorrelation];
    end
    save('FinalNeuronParameters.mat', 'NeuronParameters', 'GLM_Full_NeuronParameters', 'GLM_Partial_NeuronParameters', 'LN_NeuronParameters');
    plotResults();
end