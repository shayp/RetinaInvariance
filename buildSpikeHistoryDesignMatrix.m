function [trainspikeHistoryDesignMatrix,testspikeHistoryDesignMatrix] = buildSpikeHistoryDesignMatrix(numOfBaseVectors, numOfCoupledNeurons,...
    baseVectors, wantedLengthTrain, wantedLengthTest, spsTrain, spsTest, coupledTrain, coupledTest)

% Define the design matrix to be as the size of                                  
% coupled neurons + 1 (predicted neron)
trainspikeHistoryDesignMatrix = zeros(numOfBaseVectors * (numOfCoupledNeurons + 1), wantedLengthTrain);
testspikeHistoryDesignMatrix = zeros(numOfBaseVectors * (numOfCoupledNeurons + 1), wantedLengthTest);

% calculate the base vectors of the neuron we want to prdict
for i = 1:numOfBaseVectors  
    % make convolution of the raw stimulus(fine resolution) with a base
    % vector
    %trainspikeHistoryDesignMatrix(i,:) = conv(spsTrain,(baseVectors(:,i)'), 'same');   
    %testspikeHistoryDesignMatrix(i,:) = conv(spsTest,(baseVectors(:,i)'), 'same');
end

% calculate the  base vectors of the coupling neurons
for k = 1:numOfCoupledNeurons
    trainVector = coupledTrain(k,:);
    testVector = coupledTest(k,:);
    for i = 1:numOfBaseVectors
        
        % make convolution of the raw stimulus(fine resolution) with a base
        % vector
        trainspikeHistoryDesignMatrix(k * numOfBaseVectors + i,:) = conv(trainVector, (baseVectors(:,i)'), 'same');
        testspikeHistoryDesignMatrix(k * numOfBaseVectors + i,:) = conv(testVector, (baseVectors(:,i)'), 'same');
    end 
end

end