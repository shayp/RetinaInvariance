function [stimlusFilters, fineStimulusFilters] = getStimuliusFilterForAllCells(spikes, stimulusDesignMatrix, stimulusFilterLength)
load('globalParams');    
numOfCells = length(spikes);
    fineStimulusFilterSize = size(stimulusDesignMatrix,2);
    stimlusFilters = zeros(stimulusFilterLength, numOfCells);
    fineStimulusFilters = zeros(fineStimulusFilterSize,numOfCells);
    fineTimeScale = linspace(- deltaT * fineStimulusFilterSize,0, fineStimulusFilterSize);
    coarseTimeScale = linspace(-deltaT * fineStimulusFilterSize, 0, stimulusFilterLength);

    for i = 1:numOfCells
        tmpStimulusFilter = calculateSTA(stimulusDesignMatrix, spikes(i).data);
        fineStimulusFilters(:,i) = tmpStimulusFilter;
        stimlusFilters(:,i) = interp1(fineTimeScale, tmpStimulusFilter, coarseTimeScale, 'spline');
    end
end