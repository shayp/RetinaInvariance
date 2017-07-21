function [spikeHistoryDesignMatrix] = buildSpikeHistoryDesignMatrix(numOfBaseVectors,...
    baseVectors,spikeTrain)

    lengthOfSpikesTrain = length(spikeTrain);
    % Define the design matrix
    spikeHistoryDesignMatrix = zeros(numOfBaseVectors, lengthOfSpikesTrain);

    % calculate the  base vectors
    for i = 1:numOfBaseVectors
        % make convolution of the raw stimulus(fine resolution) with a base vector
        spikeHistoryDesignMatrix(i,:) = conv(spikeTrain, (baseVectors(:,i)'), 'same');
    end
end