function [spikes] = changeSpikeRsolution(spikeTimes, lastIndex, wantedSampFactor)
    numOfNeurons = length(spikeTimes);
    tempVector = 1:lastIndex;
    strech = ceil(lastIndex / wantedSampFactor);
    for i = 1:numOfNeurons
        spikesInSampleSize =  double(ismember(tempVector, spikeTimes(i).sp));
        spikes(i).data = zeros(strech, 1);
        for j = 1: strech  -1
            spikes(i).data(j) = sum(spikesInSampleSize((j - 1) * wantedSampFactor + 1: j * wantedSampFactor));
        end
    end
end