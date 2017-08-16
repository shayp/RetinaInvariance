function spikeHistoryData = getSpikeHistoryDataForNeurons(spikes, numOfBaseVectors, postSpikeHistory)
    numOfNeurons = length(spikes);
    
    for i = 1:numOfNeurons
        spikeHistoryData(i).spikeHistoryDesignMatrix =  buildSpikeHistoryDesignMatrix(numOfBaseVectors,...
                                                        postSpikeHistory(i).baseVectors,spikes(i).data);
    end
end