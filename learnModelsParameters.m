function neuronParameters = learnModelsParameters(neuronsInNetwork)
    load('stimulusDesignMatrix');
    load('globalParams'); 
    load('stimlusFilters');
    load('spikes');
    load('postSpikeBaseVectors');
    
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

     % get train and test spike design matrix
    [trainSpikeHistoryDesignMatrix, testSpikeHistoryDesignMatrix] = splitSpikeHistoryDesignMatrixForLearning(spikeTrainLength, neuronsInNetwork);
     
    % Get network STA filters
    realStimulusFilters = stimlusFilters(:, neuronsInNetwork);
    
    % Run for each neuron in the network
    for i = 1:numOfNeurons
        
        % Split train and test spike train
        trainSpikesTrain = spikes(neuronsInNetwork(i)).data(1:spikeTrainLength);
        testSpikesTrain = spikes(neuronsInNetwork(i)).data(spikeTrainLength + 1:end);
        
        % Get current neuron STA filter
        initStimulusFilter = realStimulusFilters(:,i);
        
        % Learn models for current neuron
        [result_GLM_Full, result_GLM_Partial, result_LN] = ...
        runLearningModels(trainStimulusDesignMatrix, testStimulusDesignMatrix,...
        trainSpikeHistoryDesignMatrix, testSpikeHistoryDesignMatrix,...
        interpMatrixTrain, interpMatrixTest, trainSpikesTrain, testSpikesTrain, stimulusFilterParamsSize,...
        binsInSecond,initStimulusFilter, numOfNeurons, numOfBaseVectors, postSpikeBaseVectors)
    
        neuronParameters(i).neuronNumber = neuronsInNetwork(i);
        neuronParameters(i).stimulusFilter = initStimulusFilter;
        neuronParameters(i).fullGLMParams = result_GLM_Full;
        neuronParameters(i).partialGLMParams = result_GLM_Partial;
        neuronParameters(i).lnOptParams = result_LN;      
    end
    
end