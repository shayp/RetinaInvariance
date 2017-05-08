function spikeHistoryDesignMatrix = buildSpikeHistoryDesignMatrix(numOfBaseVectors, numOfCoupledNeurons, baseVectors,...
    wantedLength, rawSpikesVector, wantedSampFactor, coupledVectors)

% Define the design matrix to be as the size of                                  
% coupled neurons + 1 (predicted neron)
spikeHistoryDesignMatrix = zeros(numOfBaseVectors * (numOfCoupledNeurons + 1), wantedLength);

% calculate the base vectors of the neuron we want to prdict
for i = 1:numOfBaseVectors  
    % make convolution of the raw stimulus(fine resolution) with a base
    % vector

    convVector = conv(double(rawSpikesVector), double(baseVectors(:,i))', 'same');
    
    interpSpikes = zeros(1,wantedLength);
    
    % Squeeze the result to wanted resolution
    for j = 1:wantedLength - 1
        interpSpikes(j) = sum(convVector((j - 1) * wantedSampFactor + 1: j * wantedSampFactor));
    end
    
    interpSpikes(wantedLength) = sum(convVector(j * wantedSampFactor:end));
    
    % save the squeezed convlution vector
    spikeHistoryDesignMatrix(i,:) = interpSpikes / wantedSampFactor;
end

% calculate the  base vectors of the coupling neurons
for k = 1:numOfCoupledNeurons
    rawVector = coupledVectors(k,:);
    for i = 1:numOfBaseVectors
        % make convolution of the raw stimulus(fine resolution) with a base
        % vector
        convVector = conv(double(rawVector), double(baseVectors(:,i))', 'same');
        interpSpikes = zeros(1,wantedLength);
        
        % Squeeze the result to wanted resolution
        for j = 1:wantedLength - 1
            interpSpikes(j) = sum(convVector((j - 1) * wantedSampFactor + 1: j * wantedSampFactor));
        end
        interpSpikes(wantedLength) = sum(convVector(j * wantedSampFactor:end));
        % save the squeezed convlution vector
        %spikeHistoryDesignMatrix(k * numOfBaseVectors + i,:) = interpSpikes / wantedSampFactor;
        spikeHistoryDesignMatrix(k * numOfBaseVectors + i,:) = zeros(1,wantedLength);
    end 
end
end