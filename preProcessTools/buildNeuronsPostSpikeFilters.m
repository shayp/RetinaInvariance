function postSpikeHistory = buildNeuronsPostSpikeFilters(dt, lastPeak, b, numOfVectors, absoluteRefractoryArray, ISIPeakArr)
    numOfNeurons = length(absoluteRefractoryArray);
    for i = 1:numOfNeurons
        firstPeak = absoluteRefractoryArray(i) + dt;
        peaks = [firstPeak lastPeak];
        [~, ~, baseVectors] =  buildBaseVectorsForPostSpikeAndCoupling(numOfVectors,dt,peaks, b, absoluteRefractoryArray(i));
        %lastIndex = min(75, size(baseVectors,1));
        postSpikeHistory(i).baseVectors = baseVectors;
    end
end