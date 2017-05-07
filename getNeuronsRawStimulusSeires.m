function neuronRawStimulus = getNeuronsRawStimulusSeires(numOfneurons, neurons, lastStimulus)
    neuronRawStimulus = zeros(numOfneurons,lastStimulus);
    for i = 1:numOfneurons
        neuronRawStimulus(i,:) =  ismember(1:lastStimulus, neurons(i).spikeTimes);
    end
end