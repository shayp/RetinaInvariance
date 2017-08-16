function [refreactoryPeriodArr, ISIPeakArr] =  getRefractoryPeriodForNeurons(spikes)
numOfNeurons = length(spikes);
refreactoryPeriodArr = zeros(numOfNeurons, 1);
ISIPeakArr = zeros(numOfNeurons, 1);

for i =1:numOfNeurons
    spikeTrain = spikes(i).data;
    ISI =  diff(find(spikeTrain));
    maxISI = max(ISI);
    ISIPr = zeros(maxISI, 1);
    for j = 1:maxISI
        ISIPr(j) = sum(ISI == j);
    end
    ISIPr = ISIPr / sum(ISIPr);
    [maxDiff, wantedIndex] = max(diff(ISIPr));
    if wantedIndex > 10
        wantedIndex = 10;
    end
    refreactoryPeriodArr(i) = max(2,wantedIndex - 1);
    
    [~, ISIPeakArr(i)] = max(ISIPr);
end

end