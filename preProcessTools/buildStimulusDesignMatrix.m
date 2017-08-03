function stimulusDesignMatrix = buildStimulusDesignMatrix(filterSizeBeforeSpike, stimulus)
    
    stimulusDesignMatrix = [hankel([zeros(filterSizeBeforeSpike-1,1);stimulus(1:end-filterSizeBeforeSpike+1)], ... 
        stimulus(end-filterSizeBeforeSpike+1:end))];
end