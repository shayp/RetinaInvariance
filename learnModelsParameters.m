function neuronParameters = learnModelsParameters(neuronsInNetwork)
    load('stimulusDesignMatrix');
    load('globalParams'); 
    load('stimlusFilters');
    load('spikes');
    load('postSpikeHistory');
    load('couplingBaseVectors');
    load('spikes');
    load('stimlusFilters');
    load('stimulus');
    numOfNeurons = length(neuronsInNetwork);
    
    % get train and test design matrix
    [trainStimulusDesignMatrix, testStimulusDesignMatrix] = splitStimulusDesignMatrix(stimulusDesignMatrix, trainFrac);
    stimulusTrainLength = size(trainStimulusDesignMatrix , 1);
    stimulusTestLength = size(testStimulusDesignMatrix , 1);
    spikeTrainLength = stimulusTrainLength * ceil(stimulusWantedSampleFactor / spikesWantedSampFactor);    
    spikeTestLength = stimulusTestLength * ceil(stimulusWantedSampleFactor / spikesWantedSampFactor);
    
    % Build interpolation matrix for changing resolution
    interpMatrixTrain = kron(speye(stimulusTrainLength),ones(ceil(stimulusWantedSampleFactor / spikesWantedSampFactor),1));
    interpMatrixTest = kron(speye(stimulusTestLength),ones(ceil(stimulusWantedSampleFactor / spikesWantedSampFactor),1));
     
    % Get network STA filters
    networkStimulusFilters = stimlusFilters(:, neuronsInNetwork);
    
    % Run for each neuron in the network
    for i = 1:numOfNeurons
        
        coupledNeuronnsIndexes = find(neuronsInNetwork ~= neuronsInNetwork(i));
        coupledNeurons = neuronsInNetwork(coupledNeuronnsIndexes);
        
        % get train and test spike design matrix
        [trainSpikeHistoryDesignMatrix, testSpikeHistoryDesignMatrix] = splitSpikeHistoryDesignMatrixForLearning(spikeTrainLength, neuronsInNetwork(i), coupledNeurons);
        
        % Split train and test spike train
        trainSpikesTrain = spikes(neuronsInNetwork(i)).data(1:spikeTrainLength);
        testSpikesTrain = spikes(neuronsInNetwork(i)).data(spikeTrainLength + 1:end);
        
        % Get current neuron STA filter
        initStimulusFilter = networkStimulusFilters(:,i);
        
        % Learn optimization models for current neuron
        [result_GLM_Full, result_LN] = ...
        runLearningModels(trainStimulusDesignMatrix, testStimulusDesignMatrix,...
        trainSpikeHistoryDesignMatrix, testSpikeHistoryDesignMatrix,...
        interpMatrixTrain, interpMatrixTest, trainSpikesTrain, testSpikesTrain, stimulusFilterParamsSize,...
        binsInSecond,initStimulusFilter, length(coupledNeurons), numOfBaseVectorsHistory, numOfBaseVectorsCoupling, postSpikeHistory(neuronsInNetwork(i)).baseVectors, couplingBaseVectors, stimulusFilterSizeForSimulation);

    % LN Busgang learning
        [expFunction, xData,yData] = runLNBusgang(fineStimulusFilters(:,neuronsInNetwork(i)), fineStimulus,  spikes(neuronsInNetwork(i)).data, windowSizeForFiringRate, 20);
    
        neuronParameters(i).neuronNumber = neuronsInNetwork(i);
        neuronParameters(i).coupledNeurons = coupledNeurons;
        neuronParameters(i).stimulusFilter = fineStimulusFilters(:,neuronsInNetwork(i));
        neuronParameters(i).fullGLMParams = result_GLM_Full;
        neuronParameters(i).lnOptParams = result_LN;
        neuronParameters(i).lnBusgang.expFunction = expFunction;
        neuronParameters(i).lnBusgang.xData = xData;
        neuronParameters(i).lnBusgang.yData = yData;
    end
    

end