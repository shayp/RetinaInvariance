% changeSpikeRsolution 
% We get the spike times of all neurons in the current expirment and bin
% the spike train to the new bin size
function [spikes] = changeSpikeResolution(spikeTimes, lastIndex, wantedSampFactor)

    % Get the number of neurons that have been used in the expiriment
    numOfNeurons = length(spikeTimes);
    
    tempVector = 1:lastIndex;
    
    % Calculate the new size for the spike train
    strech = ceil(lastIndex / wantedSampFactor);
    
    % Run for each neuron
    for i = 1:numOfNeurons
        % Get the fine temporal resolution of the spike train
        spikesInSampleSize =  double(ismember(tempVector, spikeTimes(i).sp));
        fineSpikeIndexes = find(spikesInSampleSize);
        spikes(i).data = zeros(strech, 1);
        
%         % Bin the spike train to more coarse resolution
%         for j = 1: strech  -1
%             spikes(i).data(j) = sum(spikesInSampleSize((j - 1) * wantedSampFactor + 1: j * wantedSampFactor));
%         end
        spikedIndexes = floor(fineSpikeIndexes / wantedSampFactor);
        spikes(i).data(spikedIndexes) = ones(length(spikedIndexes), 1);
    end
end