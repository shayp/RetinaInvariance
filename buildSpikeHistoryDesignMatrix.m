function [trainspikeHistoryDesignMatrix,testspikeHistoryDesignMatrix] = buildSpikeHistoryDesignMatrix(numOfBaseVectors, numOfCoupledNeurons,...
    baseVectors, wantedLengthTrain, wantedLengthTest, coupledTrain, coupledTest)

% Define the design matrix to be as the size of                                  
% coupled neurons + 1 (predicted neron)
trainspikeHistoryDesignMatrix = zeros(numOfBaseVectors * numOfCoupledNeurons, wantedLengthTrain);
testspikeHistoryDesignMatrix = zeros(numOfBaseVectors * numOfCoupledNeurons, wantedLengthTest);

% calculate the  base vectors of the coupling neurons
for k = 1:numOfCoupledNeurons
    trainVector = coupledTrain(k,:);
    testVector = coupledTest(k,:);
    for i = 1:numOfBaseVectors
        
        % make convolution of the raw stimulus(fine resolution) with a base
        % vector
        trainspikeHistoryDesignMatrix((k - 1) * numOfBaseVectors + i,:) = conv(trainVector, (baseVectors(:,i)'), 'same');
        testspikeHistoryDesignMatrix((k - 1) * numOfBaseVectors + i,:) = conv(testVector, (baseVectors(:,i)'), 'same');
    end 
end

end