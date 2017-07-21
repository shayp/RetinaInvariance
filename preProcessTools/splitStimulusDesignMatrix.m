function [trainStimulusDesignMatrix, testStimulusDesignMatrix] = splitStimulusDesignMatrix(stimulusDesignMatrix, trainFrac)
    stimulusLength = size(stimulusDesignMatrix, 1);
    trainLength = ceil(stimulusLength * trainFrac);
    trainStimulusDesignMatrix = stimulusDesignMatrix(1:trainLength,:);
    testStimulusDesignMatrix = stimulusDesignMatrix(trainLength + 1:end,:);
end