function [trainspikeHistoryDesignMatrix,testspikeHistoryDesignMatrix] = buildSpikeHistoryDesignMatrix(numOfBaseVectors, numOfCoupledNeurons,...
    baseVectors, wantedLengthTrain, wantedLengthTest, rawSpikesVector, wantedSampFactor, coupledVectors)

nTrainLength = wantedSampFactor * wantedLengthTrain;

% Define the design matrix to be as the size of                                  
% coupled neurons + 1 (predicted neron)
trainspikeHistoryDesignMatrix = zeros(numOfBaseVectors * (numOfCoupledNeurons + 1), wantedLengthTrain);
testspikeHistoryDesignMatrix = zeros(numOfBaseVectors * (numOfCoupledNeurons + 1), wantedLengthTest);

% calculate the base vectors of the neuron we want to prdict
for i = 1:numOfBaseVectors  
    % make convolution of the raw stimulus(fine resolution) with a base
    % vector

    convVector = conv(rawSpikesVector, fliplr(baseVectors(:,i)'), 'same');
    trainSummedSpikes = ones(1,wantedLengthTrain);
    testSummedSpikes = ones(1,wantedLengthTest);
    
    % Squeeze the result to wanted resolution
    for j = 1:wantedLengthTrain - 1
        trainSummedSpikes(j) = sum(convVector((j - 1) * wantedSampFactor + 1: j * wantedSampFactor));
    end
    for j = 1:wantedLengthTest - 1
        testSummedSpikes(j) = sum(convVector(nTrainLength + (j - 1) * wantedSampFactor + 1:nTrainLength + j * wantedSampFactor));
    end

    % save the squeezed convlution vector
    trainspikeHistoryDesignMatrix(i,:) = trainSummedSpikes;
    testspikeHistoryDesignMatrix(i,:) = testSummedSpikes;

end

% calculate the  base vectors of the coupling neurons
for k = 1:numOfCoupledNeurons
    rawVector = coupledVectors(k,:);
    for i = 1:numOfBaseVectors
        
        % make convolution of the raw stimulus(fine resolution) with a base
        % vector
        convVector = conv(rawVector, fliplr(baseVectors(:,i)'), 'same');
        trainSummedSpikes = zeros(1,wantedLengthTrain);
        testSummedSpikes = zeros(1,wantedLengthTest);        
        % Squeeze the result to wanted resolution
        for j = 1:wantedLengthTrain - 1
            trainSummedSpikes(j) = sum(convVector((j - 1) * wantedSampFactor + 1: j * wantedSampFactor));
        end
        for j = 1:wantedLengthTest - 1
            testSummedSpikes(j) = sum(convVector(nTrainLength + (j - 1) * wantedSampFactor + 1:nTrainLength + j * wantedSampFactor));
        end

        % save the squeezed convlution vector
        trainspikeHistoryDesignMatrix(k * numOfBaseVectors + i,:) = trainSummedSpikes;
        testspikeHistoryDesignMatrix(k * numOfBaseVectors + i,:) = testSummedSpikes;
    end 
end

end