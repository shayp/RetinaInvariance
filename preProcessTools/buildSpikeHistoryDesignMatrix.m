function [spikeHistoryDesignMatrix] = buildSpikeHistoryDesignMatrix(numOfBaseVectors,...
    baseVectors,spikeTrain)

    lengthOfSpikesTrain = length(spikeTrain);
    % Define the design matrix
    spikeHistoryDesignMatrix = zeros(numOfBaseVectors, lengthOfSpikesTrain);

    % calculate the  base vectors

        [lengthOFBaseVectors,~] = size(baseVectors);

        % Do convolution and remove extra bins
        Y = conv2(spikeTrain,baseVectors,'full');
        Y = [zeros(1,numOfBaseVectors); Y(1:end - lengthOFBaseVectors,:)];
        spikeHistoryDesignMatrix = Y;
end