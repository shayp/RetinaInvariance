function binnedData = binData(Data,scailingFactor)
lengthOfData = length(Data);
lengthOfBinnedData = ceil(lengthOfData / scailingFactor);
binnedData = zeros(lengthOfBinnedData,1);
    for i = 1:lengthOfBinnedData - 1
        binnedData(i) = sum(Data((i - 1) * scailingFactor + 1:i * scailingFactor));
    end
    binnedData(end) = sum(Data(i * scailingFactor:end));
    binnedData  = binnedData / scailingFactor;
end