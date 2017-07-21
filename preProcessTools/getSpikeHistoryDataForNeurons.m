function spikeHistoryData = getSpikeHistoryDataForNeurons(spikes, numOfBaseVectors, baseVectors)

    numOfNeurons = length(spikes);

    for i = 1:numOfNeurons
        spikeHistoryData(i).spikeHistoryDesignMatrix =  buildSpikeHistoryDesignMatrix(numOfBaseVectors,...
                                                        baseVectors,spikes(i).data);
    end
end