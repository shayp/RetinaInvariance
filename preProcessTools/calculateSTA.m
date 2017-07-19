function STA = calculateSTA(StimulusDesignMAtrix, cellSpikesVector)
    numOfSpikes = sum(cellSpikesVector);
    STA = (StimulusDesignMAtrix' * cellSpikesVector) ./ numOfSpikes;
end