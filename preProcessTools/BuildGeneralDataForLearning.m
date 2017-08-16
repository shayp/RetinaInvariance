function [spikes, stimulus, stimulusDesignMatrix, couplingBaseVectors, spikeHistoryData, couplingData, refreactoryPeriodArr] = BuildGeneralDataForLearning(Stim, stimtimes, SpTimes,stimulusFilterParamsSize,...
    spikesWantedSampFactor, stimulusWantedSampleFactor)
load('globalParams'); 

stimulus = changeStimulusResolution(Stim(1:end - 1),stimtimes(1:end -1), stimulusWantedSampleFactor);

lastIndex = length(stimulus) * stimulusWantedSampleFactor;

stimulusDesignMatrix = buildStimulusDesignMatrix(stimulusFilterParamsSize, stimulus);
save('stimulusDesignMatrix', 'stimulusDesignMatrix');

spikes = changeSpikeResolution(SpTimes, lastIndex, spikesWantedSampFactor);
spikesLength = length(spikes(1).data);
save('spikes', 'spikes');

[refreactoryPeriodArr,ISIPeakArr] =  getRefractoryPeriodForNeurons(spikes);
refreactoryPeriodArr = refreactoryPeriodArr * dt;
ISIPeakArr = ISIPeakArr * dt;
save('ISIPeakArr', 'ISIPeakArr');
save('refreactoryPeriodArr', 'refreactoryPeriodArr');
% Get the fine stimulus resolution(Just for calculating the STA)
fineStimulus = changeStimulusResolution(Stim,stimtimes, spikesWantedSampFactor);

fineStimulus = fineStimulus(1:spikesLength);
save('stimulus', 'stimulus', 'fineStimulus');

% Calculate the fine stimulus DesignMatrix(For STA only)
fineStimulusDesignMatrix = buildStimulusDesignMatrix(stimulusFilterParamsSize * (stimulusWantedSampleFactor / spikesWantedSampFactor), fineStimulus);

[stimlusFilters, fineStimulusFilters] = getStimuliusFilterForAllCells(spikes, fineStimulusDesignMatrix, stimulusFilterParamsSize);
save('stimlusFilters', 'stimlusFilters', 'fineStimulusFilters');

postSpikeHistory = buildNeuronsPostSpikeFilters(dt, lastPeakHistory, b, numOfBaseVectorsHistory, refreactoryPeriodArr, ISIPeakArr);
save('postSpikeHistory','postSpikeHistory');
% build coupling BaseVectors
[~, ~, couplingBaseVectors] = buildBaseVectorsForPostSpikeAndCoupling(numOfBaseVectorsCoupling ,dt, [dt lastPeakCoupling], b, 0);
 save('couplingBaseVectors', 'couplingBaseVectors');

% 
% figure();
% plot(couplingBaseVectors);drawnow;

spikeHistoryData = getSpikeHistoryDataForNeurons(spikes, numOfBaseVectorsHistory, postSpikeHistory);
save('spikeHistoryData', 'spikeHistoryData', '-v7.3');

couplingData = getCouplingDataForNeurons(spikes, numOfBaseVectorsCoupling, couplingBaseVectors);
save('couplingData', 'couplingData', '-v7.3');

end