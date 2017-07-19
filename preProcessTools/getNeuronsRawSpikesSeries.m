function [neuronRawSpikes, neuronsSclaedSpikes] = getNeuronsRawSpikesSeries(numOfneurons, neurons, lastStimulus, wantedSampFactor, scaledSize)
    neuronRawSpikes = zeros(numOfneurons,lastStimulus);
    neuronsSclaedSpikes = zeros(numOfneurons,scaledSize);
    for i = 1:numOfneurons
        neuronRawSpikes(i,:) =  double(ismember(1:lastStimulus, neurons(i).spikeTimes));
        for j = 1:scaledSize - 1
            neuronsSclaedSpikes(i,j) = sum(neuronRawSpikes(i,(j - 1) * wantedSampFactor + 1: j * wantedSampFactor));
        end
    end
end