function postSpikeHistory = buildNeuronsPostSpikeFilters(dt, lastPeak, b, numOfVectors, absoluteRefractoryArray, ISIPeakArr)
    numOfNeurons = length(absoluteRefractoryArray);
    for i = 1:numOfNeurons
        firstPeak = max(absoluteRefractoryArray(i) + dt, ISIPeakArr(i));
        peaks = [firstPeak lastPeak];
        [~, ~, postSpikeHistory(i).baseVectors] =  buildBaseVectorsForPostSpikeAndCoupling(numOfVectors,dt,peaks, b, absoluteRefractoryArray(i));
    end
end