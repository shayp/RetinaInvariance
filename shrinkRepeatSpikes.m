function scaledRepSpikes = shrinkRepeatSpikes(stimTimes, spikeTimes, wantedSampFactor)
    vectorLen = ceil((stimTimes(2) - stimTimes(1)) /  wantedSampFactor);
    numOfReturns = length(stimTimes);
    scaledRepSpikes = zeros(numOfReturns,vectorLen);
    vectorLen
    for i = 1:numOfReturns - 1
        firstIndex = stimTimes(i);
        lastIndex = stimTimes(i + 1) - 1;
        currentReturnSpike = double(ismember(firstIndex:lastIndex, spikeTimes));
        currentReturnSpike
        for j = 1:vectorLen - 1
            scaledRepSpikes(i,j) = sum(currentReturnSpike((j - 1) * wantedSampFactor + 1: j * wantedSampFactor));
        end
    end
end