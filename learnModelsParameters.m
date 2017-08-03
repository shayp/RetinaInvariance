function neuronParameters = learnModelsParameters(neuronsInNetwork)
    load('stimulusDesignMatrix');
    load('globalParams'); 
    load('stimlusFilters');
    load('spikes');
    load('postSpikeBaseVectors');
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
        
        % Learn optimization models for current neuron
        [result_GLM_Full, result_LN] = ...
        runLearningModels(trainStimulusDesignMatrix, testStimulusDesignMatrix,...
        trainSpikeHistoryDesignMatrix, testSpikeHistoryDesignMatrix,...
        interpMatrixTrain, interpMatrixTest, trainSpikesTrain, testSpikesTrain, stimulusFilterParamsSize,...
        binsInSecond,initStimulusFilter, numOfNeurons, numOfBaseVectors, postSpikeBaseVectors, stimulusFilterSizeForSimulation);
    
        % LN Busgang learning
        [expFunction, xData,yData] = runLNBusgang(fineStimulusFilters(:,neuronsInNetwork(i)), fineStimulus,  spikes(neuronsInNetwork(i)).data, windowSizeForFiringRate, 20);
    
        neuronParameters(i).neuronNumber = neuronsInNetwork(i);
        neuronParameters(i).stimulusFilter = fineStimulusFilters(:,neuronsInNetwork(i));
        neuronParameters(i).fullGLMParams = result_GLM_Full;
        %neuronParameters(i).partialGLMParams = result_GLM_Partial;
        neuronParameters(i).lnOptParams = result_LN;
        neuronParameters(i).lnBusgang.expFunction = expFunction;
        neuronParameters(i).lnBusgang.xData = xData;
        neuronParameters(i).lnBusgang.yData = yData;
    end
    
% timeSeries = linspace(-stimulusFilterSizeForSimulation * deltaT, 0, stimulusFilterSizeForSimulation);
% figure();
% plot(timeSeries, neuronParameters(1).stimulusFilter,...
%      timeSeries, neuronParameters(1).lnOptParams.StimulusFilter,...
%      timeSeries, neuronParameters(1).fullGLMParams.StimulusFilter,...
%      timeSeries, neuronParameters(1).partialGLMParams.StimulusFilter);
% legend('STA', 'LN', 'GLM Full', 'GLM Partial');
% xlabel('Time before spike(s)');
% ylabel('intensity');
% title('Choosen stimulus filter(after resize):');
% drawnow;
end