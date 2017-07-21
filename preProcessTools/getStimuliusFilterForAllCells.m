function [stimlusFilters, fineStimulusFilters] = getStimuliusFilterForAllCells(spikes, stimulusDesignMatrix, stimulusFilterLength)
    numOfCells = length(spikes);
    fineStimulusFilterSize = size(stimulusDesignMatrix,2);
    stimlusFilters = zeros(stimulusFilterLength, numOfCells);
    fineStimulusFilters = zeros(fineStimulusFilterSize,numOfCells);
    for i = 1:numOfCells
        tmpStimulusFilter = calculateSTA(stimulusDesignMatrix, spikes(i).data);
        fineStimulusFilters(:,i) = tmpStimulusFilter;
        stimlusFilters(:,i) = imresize(tmpStimulusFilter, [stimulusFilterLength 1]);
    end
end