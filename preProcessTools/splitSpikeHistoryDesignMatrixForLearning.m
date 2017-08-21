function [trainSpikeHistoryDesignMatrix, testSpikeHistoryDesignMatrix] = splitSpikeHistoryDesignMatrixForLearning(trainLength,selectedNeuron, coupledNeurons)
load('spikeHistoryData');
load('couplingData');
load('globalParams');

numOfCoupledNeurons =  length(coupledNeurons);
tempDesignMatrix = spikeHistoryData(1).spikeHistoryDesignMatrix;
spikeLength = size(tempDesignMatrix,1);
testLength = spikeLength - trainLength;

trainSpikeHistoryDesignMatrix = zeros(trainLength, numOfBaseVectorsHistory + numOfBaseVectorsCoupling * numOfCoupledNeurons);
testSpikeHistoryDesignMatrix = zeros(testLength, numOfBaseVectorsHistory + numOfBaseVectorsCoupling * numOfCoupledNeurons);

spikeHistoryDesignMatrix = spikeHistoryData(selectedNeuron).spikeHistoryDesignMatrix;
trainSpikeHistoryDesignMatrix(:, 1:numOfBaseVectorsHistory) = spikeHistoryDesignMatrix(1:trainLength,:);
testSpikeHistoryDesignMatrix(:, 1:numOfBaseVectorsHistory) = spikeHistoryDesignMatrix(trainLength + 1:end,:);
    
for i = 1:numOfCoupledNeurons
    currentDesignMatrix = couplingData(coupledNeurons(i)).couplingDesignMatrix;
    trainSpikeHistoryDesignMatrix(:,numOfBaseVectorsHistory + (i - 1) * numOfBaseVectorsCoupling + 1:numOfBaseVectorsHistory + i * numOfBaseVectorsCoupling) = ...
        currentDesignMatrix(1:trainLength,:);
    testSpikeHistoryDesignMatrix(:,numOfBaseVectorsHistory + (i - 1) * numOfBaseVectorsCoupling + 1:numOfBaseVectorsHistory + i * numOfBaseVectorsCoupling) = ...
        currentDesignMatrix(trainLength + 1:end,:);
end

end