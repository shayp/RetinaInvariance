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
     wantedIndexes = find(ISIPr >= 0.01);
     if isempty(wantedIndexes)
         wantedIndex = 2;
     else
        wantedIndex = wantedIndexes(1) - 1;
     end
     
    if wantedIndex > 4
        wantedIndex = 4;
    end
    refreactoryPeriodArr(i) = max(2,wantedIndex);
    
    [~, ISIPeak] = max(ISIPr);
    ISIPeakArr(i) = min(10, ISIPeak);
end

end