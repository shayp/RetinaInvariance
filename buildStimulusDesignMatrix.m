function stimulusDesignMatrix = buildStimulusDesignMatrix(filterSizeBeforeSpike, stimulus)
    stimLength = length(stimulus);
    
    stimulusDesignMatrix = [hankel([zeros(filterSizeBeforeSpike-1,1);stimulus(1:end-filterSizeBeforeSpike+1)], ... 
        stimulus(end-filterSizeBeforeSpike+1:end))];
end