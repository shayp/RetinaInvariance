function couplingData = getCouplingDataForNeurons(spikes, numOfBaseVectors, couplingBaseVectors)
    numOfNeurons = length(spikes);
    
    for i = 1:numOfNeurons
        couplingData(i).couplingDesignMatrix =  buildSpikeHistoryDesignMatrix(numOfBaseVectors,...
                                                       couplingBaseVectors,spikes(i).data);
    end
end