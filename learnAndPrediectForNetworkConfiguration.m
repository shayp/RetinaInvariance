function learnAndPrediectForNetworkConfiguration(neuronsInNetwork, scaledRepStimulus)
    load( 'SpTimes'); 
    load ('repeatStimulusTimes');    
    load ('RepStimulusExtended');
    load ('RepSpTimes'); 
    load('globalParams');
    

    numOfNeurons = length(neuronsInNetwork);
    learnedParameters = learnModelsParameters(neuronsInNetwork);
    NeuronParameters = getRepStimulusData(neuronsInNetwork, repeatStimulusTimes, RepSpTimes, spikesWantedSampFactor, RepStimulusExtended);
    
    for i = 1:numOfNeurons
        NeuronParameters(i).neuronsInNetwork = neuronsInNetwork;
        
        NeuronParameters(i).neuronIndex = learnedParameters(i).neuronNumber;
        glmFullParams(i).neuronIndex = learnedParameters(i).neuronNumber;
        glmPartialParams(i).neuronIndex = learnedParameters(i).neuronNumber;
        lnParams(i).neuronIndex = learnedParameters(i).neuronNumber;

        NeuronParameters(i).stimulusFilter = learnedParameters(i).stimulusFilter;
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
        
        lnBusgangParams(i).stimulusFilter = learnedParameters(i).stimulusFilter;
        lnBusgangParams(i).expFunction = learnedParameters(i).lnBusgang.expFunction;
        lnBusgangParams(i).xData = learnedParameters(i).lnBusgang.xData;
        lnBusgangParams(i).yData = learnedParameters(i).lnBusgang.yData;

    end

    for j = 1:numOfRepeats
        
        % Full glm simulation
        response_GLM_Full = RunSimulation(numOfNeurons, scaledRepStimulus, glmFullParams, stimulusFilterSizeForSimulation, baseVectorLength, deltaT, 1);
        responseGLM_GLM_Partial = RunSimulation(numOfNeurons, scaledRepStimulus, glmPartialParams, stimulusFilterSizeForSimulation, baseVectorLength, deltaT , 1);
        response_LN = RunSimulation(numOfNeurons, scaledRepStimulus, lnParams, stimulusFilterSizeForSimulation, baseVectorLength, deltaT, 0);
       
        for i = 1:numOfNeurons
            NeuronParameters(i).GLMFullSimulation(j,:) =  response_GLM_Full(i,:);
            NeuronParameters(i).GLMPartialSimulation(j,:) =  responseGLM_GLM_Partial(i,:);
            NeuronParameters(i).LNSimulation(j,:) =  response_LN(i,:);
        end
    end

    numbrOfBins = 10;
    for i = 1:numOfNeurons
        NeuronParameters(i).lnBusgangFiringRate = estimateLNBugangFiringRate(scaledRepStimulus, lnBusgangParams(i).stimulusFilter, lnBusgangParams(i).expFunction, windowSizeForFiringRate, deltaT);
        
        [glmFullSpikeRate, glmFullCorrelation] = CalculateCorrelatedSpikeRate(numOfRepeats, NeuronParameters(i).scaledRepSpikes, NeuronParameters(i).GLMFullSimulation, windowSizeForFiringRate);
        [glmPartialSpikeRate, glmPartialCorrelation] = CalculateCorrelatedSpikeRate(numOfRepeats, NeuronParameters(i).scaledRepSpikes, NeuronParameters(i).GLMPartialSimulation, windowSizeForFiringRate);
        [lnSpikeRate,lnCorrelation] = CalculateCorrelatedSpikeRate(numOfRepeats, NeuronParameters(i).scaledRepSpikes, NeuronParameters(i).LNSimulation, windowSizeForFiringRate);
        [partSpikeRate, autoCorrelation] = CalculateCorrelatedSpikeRate(5, NeuronParameters(i).scaledRepSpikes(1:5,:),NeuronParameters(i).scaledRepSpikes(6:end,:), windowSizeForFiringRate);
        [lnBusgangSpikeRate, lnBusgangcorrelation] = CalculateCorrelatedSpikeRateBusgang(NeuronParameters(i).scaledRepSpikes, NeuronParameters(i).lnBusgangFiringRate, windowSizeForFiringRate);
       
        NeuronParameters(i).realSpikeRate = glmFullSpikeRate(1,:);
        NeuronParameters(i).glmFullSimulatedSpikeRate = glmFullSpikeRate(2,:);
        NeuronParameters(i).glmPartialSimulatedSpikeRate = glmPartialSpikeRate(2,:);
        NeuronParameters(i).lnSimulatedSpikeRate = lnSpikeRate(2,:);
        NeuronParameters(i).partSpikeRate = partSpikeRate(2,:);
        NeuronParameters(i).lnBusgangSpikeRate = lnBusgangSpikeRate(2,:);
        
        correaltionVector = zeros(6, numbrOfBins);
        maxRealSpike = max(NeuronParameters(i).realSpikeRate);
        firingRateSpace = linspace(0, maxRealSpike, numbrOfBins + 1);
        for bin = 1:numbrOfBins
            wantedIndexes = find(NeuronParameters(i).realSpikeRate >= firingRateSpace(bin) & NeuronParameters(i).realSpikeRate < firingRateSpace(bin + 1));
            sizeOfBin = length(wantedIndexes);
            correaltionVector(1, bin) = sum(NeuronParameters(i).realSpikeRate(wantedIndexes)) / sizeOfBin;
            correaltionVector(2, bin) = sum(NeuronParameters(i).lnSimulatedSpikeRate(wantedIndexes)) / sizeOfBin;
            correaltionVector(3, bin) = sum(NeuronParameters(i).glmPartialSimulatedSpikeRate(wantedIndexes)) / sizeOfBin;
            correaltionVector(4, bin) = sum(NeuronParameters(i).glmFullSimulatedSpikeRate(wantedIndexes)) / sizeOfBin;
            correaltionVector(5, bin) = sum(NeuronParameters(i).lnBusgangSpikeRate(wantedIndexes)) / sizeOfBin;
            correaltionVector(6, bin) = sum(NeuronParameters(i).partSpikeRate(wantedIndexes)) / sizeOfBin;
        end
        NeuronParameters(i).correaltionVector = correaltionVector;
        NeuronParameters(i).spikeRateCorrelation = [lnCorrelation glmPartialCorrelation glmFullCorrelation lnBusgangcorrelation autoCorrelation];
        lnExplained = min(lnCorrelation/ autoCorrelation * 100, 100);
        glmPartialExplained = min(glmPartialCorrelation/ autoCorrelation * 100, 100);
        glmFullExplained = min(glmFullCorrelation/ autoCorrelation * 100, 100);
        lnBusgangExplained = min(lnBusgangcorrelation/ autoCorrelation * 100, 100);

        NeuronParameters(i).perecentExplained = [lnExplained glmPartialExplained glmFullExplained lnBusgangExplained];
    end
    save('FinalNeuronParameters.mat', 'NeuronParameters', 'glmFullParams', 'glmPartialParams', 'lnParams', 'lnBusgangParams');
    plotResults(length(neuronsInNetwork));
end