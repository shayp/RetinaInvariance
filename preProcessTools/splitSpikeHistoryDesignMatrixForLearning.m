function [trainSpikeHistoryDesignMatrix, testSpikeHistoryDesignMatrix] = splitSpikeHistoryDesignMatrixForLearning(trainLength, neuronsToInclude)
load('spikeHistoryData');
numOfNeurons =  length(neuronsToInclude);
tempDesignMatrix = spikeHistoryData(1).spikeHistoryDesignMatrix;
spikeLength = size(tempDesignMatrix,2);
numOfBaseVectors = size(tempDesignMatrix,1);
testLength = spikeLength - trainLength;
trainSpikeHistoryDesignMatrix = zeros(numOfBaseVectors * numOfNeurons,trainLength);
testSpikeHistoryDesignMatrix = zeros(numOfBaseVectors * numOfNeurons,testLength);

for i = 1:numOfNeurons
    currentDesignMatrix = spikeHistoryData(neuronsToInclude(i)).spikeHistoryDesignMatrix;
    trainSpikeHistoryDesignMatrix((i - 1) * numOfBaseVectors + 1:i * numOfBaseVectors,:) = ...
        currentDesignMatrix(:,1:trainLength);
    testSpikeHistoryDesignMatrix((i - 1) * numOfBaseVectors + 1:i * numOfBaseVectors,:) = ...
        currentDesignMatrix(:,trainLength + 1:end);
end
end