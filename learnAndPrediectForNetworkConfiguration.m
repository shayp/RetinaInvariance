function learnAndPrediectForNetworkConfiguration(neuronsInNetwork, scaledRepStimulus)
    load( 'SpTimes'); 
    load ('repeatStimulusTimes');    
    load ('RepStimulusExtended');
    load ('RepSpTimes'); 
    load('globalParams');
    

    numOfNeurons = length(neuronsInNetwork);
    learnedParameters = learnModelsParameters(neuronsInNetwork);
    NeuronParameters = getRepStimulusData(neuronsInNetwork, repeatStimulusTimes, RepSpTimes, spikesWantedSampFactor);
    repStimulusDesignMatrix = buildStimulusDesignMatrix(stimulusFilterSizeForSimulation, scaledRepStimulus);
    
    for i = 1:numOfNeurons
        NeuronParameters(i).neuronsInNetwork = neuronsInNetwork;
        NeuronParameters(i).repSTA = getRepeatSTA(repStimulusDesignMatrix,NeuronParameters(i).scaledRepSpikes, stimulusFilterSizeForSimulation);
        NeuronParameters(i).neuronIndex = learnedParameters(i).neuronNumber;
        glmFullParams(i).neuronIndex = learnedParameters(i).neuronNumber;
        lnParams(i).neuronIndex = learnedParameters(i).neuronNumber;

        NeuronParameters(i).stimulusFilter = learnedParameters(i).stimulusFilter;

        glmFullParams(i).stimulusFilter = learnedParameters(i).fullGLMParams.StimulusFilter;
        lnParams(i).stimulusFilter = learnedParameters(i).lnOptParams.StimulusFilter;

        glmFullParams(i).meanFiringRate = learnedParameters(i).fullGLMParams.meanFiringRate;
        lnParams(i).meanFiringRate = learnedParameters(i).lnOptParams.meanFiringRate;

        glmFullParams(i).couplingFilters = learnedParameters(i).fullGLMParams.couplingFilters;
    
        glmFullParams(i).simulation = zeros(numOfRepeats, length(scaledRepStimulus));
        lnParams(i).simulation = zeros(numOfRepeats, length(scaledRepStimulus));
        
        lnBusgangParams(i).stimulusFilter = learnedParameters(i).stimulusFilter;
        lnBusgangParams(i).expFunction = learnedParameters(i).lnBusgang.expFunction;
        lnBusgangParams(i).xData = learnedParameters(i).lnBusgang.xData;
        lnBusgangParams(i).yData = learnedParameters(i).lnBusgang.yData;

    end

    for j = 1:numOfRepeats
        
        % Full glm simulation
        response_GLM_Full = RunSimulation(numOfNeurons, scaledRepStimulus, glmFullParams, stimulusFilterSizeForSimulation, deltaT, 1);
        response_LN = RunSimulation(numOfNeurons, scaledRepStimulus, lnParams, stimulusFilterSizeForSimulation, deltaT, 0);
       
        for i = 1:numOfNeurons
            NeuronParameters(i).GLMFullSimulation(j,:) =  response_GLM_Full(i,:);
            NeuronParameters(i).LNSimulation(j,:) =  response_LN(i,:);
        end
    end

    pathName = 'Results/Network';
    for i = 1:numOfNeurons
        NeuronParameters(i).GLMFullISI = [];
        NeuronParameters(i).LNISI = [];
        NeuronParameters(i).repISI = [];

        for j = 1:numOfRepeats
            NeuronParameters(i).GLMFullISI = [NeuronParameters(i).GLMFullISI diff(find(NeuronParameters(i).GLMFullSimulation(j,:)))];
            NeuronParameters(i).LNISI = [NeuronParameters(i).LNISI diff(find(NeuronParameters(i).LNSimulation(j,:)))];
        end
        
        for j = 1:size(NeuronParameters(i).scaledRepSpikes,1)
            NeuronParameters(i).repISI = [NeuronParameters(i).repISI diff(find(NeuronParameters(i).scaledRepSpikes(j,:)))];
        end
        pathName = strcat(pathName, ['_' num2str(NeuronParameters(i).neuronIndex)]);
    end
    

    numbrOfBins = 10;
    for i = 1:numOfNeurons
        NeuronParameters(i).lnBusgangFiringRate = estimateLNBugangFiringRate(scaledRepStimulus, lnBusgangParams(i).stimulusFilter, lnBusgangParams(i).expFunction, windowSizeForFiringRate, deltaT);
        
        [partSpikeRate, autoCorrelation, partVariance] = CalculateCorrelatedSpikeRate(5, NeuronParameters(i).scaledRepSpikes(1:5,:),NeuronParameters(i).scaledRepSpikes(6:end,:), windowSizeForFiringRate, deltaT);
        [glmFullSpikeRate, glmFullCorrelation, glmFullVariance] = CalculateCorrelatedSpikeRate(numOfRepeats, NeuronParameters(i).scaledRepSpikes, NeuronParameters(i).GLMFullSimulation, windowSizeForFiringRate, deltaT);
        [lnSpikeRate,lnCorrelation, LnVariance] = CalculateCorrelatedSpikeRate(numOfRepeats, NeuronParameters(i).scaledRepSpikes, NeuronParameters(i).LNSimulation, windowSizeForFiringRate, deltaT);
        [lnBusgangSpikeRate, lnBusgangcorrelation, LNBusgangVariance] = CalculateCorrelatedSpikeRateBusgang(NeuronParameters(i).scaledRepSpikes, NeuronParameters(i).lnBusgangFiringRate, windowSizeForFiringRate, deltaT);
       
        NeuronParameters(i).realSpikeRate = glmFullSpikeRate(1,:);
        NeuronParameters(i).glmFullSimulatedSpikeRate = glmFullSpikeRate(2,:);
        NeuronParameters(i).lnSimulatedSpikeRate = lnSpikeRate(2,:);
        NeuronParameters(i).partSpikeRate = partSpikeRate(2,:);
        NeuronParameters(i).lnBusgangSpikeRate = lnBusgangSpikeRate(2,:);
        
        
        
        correaltionVector = zeros(5, numbrOfBins);
        maxRealSpike = max(NeuronParameters(i).realSpikeRate);
        firingRateSpace = linspace(0, maxRealSpike, numbrOfBins + 1);
        for bin = 1:numbrOfBins
            wantedIndexes = find(NeuronParameters(i).realSpikeRate >= firingRateSpace(bin) & NeuronParameters(i).realSpikeRate < firingRateSpace(bin + 1));
            sizeOfBin = length(wantedIndexes);
            correaltionVector(1, bin) = sum(NeuronParameters(i).realSpikeRate(wantedIndexes)) / sizeOfBin;
            correaltionVector(2, bin) = sum(NeuronParameters(i).lnSimulatedSpikeRate(wantedIndexes)) / sizeOfBin;
            correaltionVector(3, bin) = sum(NeuronParameters(i).glmFullSimulatedSpikeRate(wantedIndexes)) / sizeOfBin;
            correaltionVector(4, bin) = sum(NeuronParameters(i).lnBusgangSpikeRate(wantedIndexes)) / sizeOfBin;
            correaltionVector(5, bin) = sum(NeuronParameters(i).partSpikeRate(wantedIndexes)) / sizeOfBin;
        end
        NeuronParameters(i).correaltionVector = correaltionVector;
        NeuronParameters(i).varianceExplained = [LnVariance glmFullVariance LNBusgangVariance partVariance]; 
        NeuronParameters(i).spikeRateCorrelation = [lnCorrelation glmFullCorrelation lnBusgangcorrelation autoCorrelation];
        lnExplained = min(lnCorrelation/ autoCorrelation * 100, 100);
        %glmPartialExplained = min(glmPartialCorrelation/ autoCorrelation * 100, 100);
        glmFullExplained = min(glmFullCorrelation/ autoCorrelation * 100, 100);
        lnBusgangExplained = min(lnBusgangcorrelation/ autoCorrelation * 100, 100);

        NeuronParameters(i).perecentExplained = [lnExplained glmFullExplained lnBusgangExplained];
    end
    pathName = strcat(pathName,'/');
    mkdir (pathName);
    save([pathName 'FinalNeuronParameters'], 'NeuronParameters', 'glmFullParams', 'lnParams', 'lnBusgangParams');
    plotResults(length(neuronsInNetwork), pathName);
end