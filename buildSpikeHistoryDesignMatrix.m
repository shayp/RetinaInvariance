function spikeHistoryDesignMatrix = buildSpikeHistoryDesignMatrix(numOfBaseVectors, baseVectors, cellSpikesVector)
    spikeHistoryDesignMatrix = zeros(numOfBaseVectors, length(cellSpikesVector));
    for i = 1:numOfBaseVectors
        spikeHistoryDesignMatrix(i,:) = conv(cellSpikesVector, baseVectors(i), 'same');
    end
end