function [spikeRate, correlation, varExplain] = CalculateCorrelatedSpikeRateBusgang(realNeuronData, simulatedNeuronData, windowSize, dt)
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
spikeRate(1,:) = spikeRate(1,:) /repeatInExp;
spikeRate = spikeRate / (windowSize * dt);
[vecCorrelation, vecLegs] = xcorr(spikeRate(1,:),spikeRate(2,:));
[~, index] = max(vecCorrelation);
Leg =  vecLegs(index)
cutOut = ceil(spikeRateLength / 30)
spikeRate(2,:) = circshift(spikeRate(2,:)',Leg)';
correlation = corr2(spikeRate(1,:),spikeRate(2,:))
correlation = corr2(spikeRate(1,cutOut:end - cutOut),spikeRate(2,cutOut:end - cutOut))

meanOfSpikesTrain =  mean(spikeRate(1,:));
totalSumOfErrors = sum((spikeRate(2,:) - spikeRate(1,:)).^2);
totalSumOfSquraes = sum((spikeRate(2,:) -meanOfSpikesTrain).^2);
varExplain = 1 -totalSumOfErrors/totalSumOfSquraes;