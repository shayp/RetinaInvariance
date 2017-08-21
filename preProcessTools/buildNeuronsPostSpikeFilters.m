function postSpikeHistory = buildNeuronsPostSpikeFilters(dt, lastPeak, b, numOfVectors, absoluteRefractoryArray, ISIPeakArr)
    numOfNeurons = length(absoluteRefractoryArray);
    for i = 1:numOfNeurons
        firstPeak = max(absoluteRefractoryArray(i) + dt, ISIPeakArr(i));
        %firstPeak = absoluteRefractoryArray(i) + dt;
        peaks = [firstPeak lastPeak];
        [~, ~, baseVectors] =  buildBaseVectorsForPostSpikeAndCoupling(numOfVectors,dt,peaks, b, absoluteRefractoryArray(i));
        postSpikeHistory(i).baseVectors = baseVectors;
    end
end