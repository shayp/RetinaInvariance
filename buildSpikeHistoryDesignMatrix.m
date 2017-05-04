function spikeHistoryDesignMatrix = buildSpikeHistoryDesignMatrix(numOfBaseVectors, baseVectors, wantedLength, rawSpikesVector, wantedSampFactor)
    spikeHistoryDesignMatrix = zeros(numOfBaseVectors, wantedLength);
    for i = 1:numOfBaseVectors  
        convVector = conv(double(rawSpikesVector), double(baseVectors(:,i))', 'same');
        interpSpikes = zeros(1,wantedLength);
        for j = 1:wantedLength  -1
            sum(convVector((j - 1) * wantedSampFactor + 1: j * wantedSampFactor));
            interpSpikes(j) = sum(convVector((j - 1) * wantedSampFactor + 1: j * wantedSampFactor));
        end
        spikeHistoryDesignMatrix(i,:) = interpSpikes;
    end
end