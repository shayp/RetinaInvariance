function stimulusDesignMatrix = buildStimulusDesignMatrix(filterSizeBeforeSpike, stimulus)
    stimLength = len(stimulus);
    % build the design matrix, training data
    stimulusDesignMatrix = [ones(stimLength,1), ... 
        hankel([zeros(filterSizeBeforeSpike-1,1);stimulus(1:end-filterSizeBeforeSpike+1)], ... 
        stimulus(end-filterSizeBeforeSpike+1:end))];
end