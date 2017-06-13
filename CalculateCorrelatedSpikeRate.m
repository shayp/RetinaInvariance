function [spikeRate, correlation] = CalculateCorrelatedSpikeRate(numOfRepeats, realNeuronData, simulatedNeuronData, windowSize)
dataLength = size(realNeuronData,2);
spikeRateLength = ceil(dataLength / windowSize);
spikeRate = zeros(2, spikeRateLength);
for j = 1:spikeRateLength - 1
    for k = 1:numOfRepeats
    spikeRate(1,j) = spikeRate(1,j) + sum(realNeuronData(k, (j-1) * windowSize + 1: j * windowSize));
    spikeRate(2,j) = spikeRate(2,j) + sum(simulatedNeuronData(k, (j-1) * windowSize + 1: j * windowSize));
    end
end
sum(spikeRate(1,:))
sum(spikeRate(2,:))
spikeRate(1,:) = circshift(spikeRate(1,:)',1)';
spikeRate = spikeRate / (windowSize * numOfRepeats);
correlation = corr2(spikeRate(1,:),spikeRate(2,:))