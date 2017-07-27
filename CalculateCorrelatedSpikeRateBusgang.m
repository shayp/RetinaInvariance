function [spikeRate, correlation] = CalculateCorrelatedSpikeRateBusgang(realNeuronData, simulatedNeuronData, windowSize)
dataLength = size(realNeuronData,2);
spikeRateLength = length(simulatedNeuronData);
spikeRate = zeros(2, spikeRateLength);
repeatInExp = size(realNeuronData, 1);
spikeRate(2,:) = simulatedNeuronData;
for j = 1:spikeRateLength
    currentChange = min(dataLength, j * windowSize);
    for i = 1:repeatInExp
       spikeRate(1,j) = spikeRate(1,j) + sum(realNeuronData(i, (j-1) * windowSize + 1: currentChange));
    end
end
spikeRate(1,:) = spikeRate(1,:) / (windowSize * repeatInExp);
spikeRate(2,:) = spikeRate(2,:) / (windowSize);
[vecCorrelation, vecLegs] = xcorr(spikeRate(1,:),spikeRate(2,:));
[~, index] = max(vecCorrelation);
Leg =  vecLegs(index)
cutOut = ceil(spikeRateLength / 30)
spikeRate(2,:) = circshift(spikeRate(2,:)',Leg)';
correlation = corr2(spikeRate(1,:),spikeRate(2,:))
correlation = corr2(spikeRate(1,cutOut:end - cutOut),spikeRate(2,cutOut:end - cutOut))
spikeRate(2,:) = spikeRate(2,:) * (max(spikeRate(1,:)) / max(spikeRate(2,:)));
