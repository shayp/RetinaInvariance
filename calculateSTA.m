function STA = calculateSTA(StimulusDesignMAtrix, cellSpikesVector)
    STA = StimulusDesignMAtrix * cellSpikes;
end