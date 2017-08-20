function [spikes, stimulus, stimulusDesignMatrix, couplingBaseVectors, spikeHistoryData, couplingData, refreactoryPeriodArr] = BuildGeneralDataForLearning(Stim, stimtimes, SpTimes,stimulusFilterParamsSize,...
    spikesWantedSampFactor, stimulusWantedSampleFactor)
load('globalParams'); 

disp('**********      changeStimulusResolution(stimulusWantedSampleFactor)      **********');
stimulus = changeStimulusResolution(Stim(1:end - 1),stimtimes(1:end -1), stimulusWantedSampleFactor);
save('stimulus','stimulus');
lastIndex = length(stimulus) * stimulusWantedSampleFactor;

disp('**********      buildStimulusDesignMatrix(coarse stimulus)      **********');
stimulusDesignMatrix = buildStimulusDesignMatrix(stimulusFilterParamsSize, stimulus);
save('stimulusDesignMatrix', 'stimulusDesignMatrix');

disp('**********      changeSpikeResolution      **********');
spikes = changeSpikeResolution(SpTimes, lastIndex, spikesWantedSampFactor);
spikesLength = length(spikes(1).data);
save('spikes', 'spikes');

disp('**********      getRefractoryPeriodForNeurons      **********');
[refreactoryPeriodArr,ISIPeakArr] =  getRefractoryPeriodForNeurons(spikes);
refreactoryPeriodArr = refreactoryPeriodArr * dt;
ISIPeakArr = ISIPeakArr * dt;
save('ISIPeakArr', 'ISIPeakArr');
save('refreactoryPeriodArr', 'refreactoryPeriodArr');

disp('**********      changeStimulusResolution(spikesWantedSampFactor)      **********');

% Get the fine stimulus resolution(Just for calculating the STA)
fineStimulus = changeStimulusResolution(Stim,stimtimes, spikesWantedSampFactor);

fineStimulus = fineStimulus(1:spikesLength);
save('stimulus', 'stimulus', 'fineStimulus');

disp('**********      buildStimulusDesignMatrix(fine stimulus)      **********');
% Calculate the fine stimulus DesignMatrix(For STA only)
fineStimulusDesignMatrix = buildStimulusDesignMatrix(stimulusFilterParamsSize * (stimulusWantedSampleFactor / spikesWantedSampFactor), fineStimulus);

disp('**********      getStimuliusFilterForAllCells      **********');
[stimlusFilters, fineStimulusFilters] = getStimuliusFilterForAllCells(spikes, fineStimulusDesignMatrix, stimulusFilterParamsSize);
save('stimlusFilters', 'stimlusFilters', 'fineStimulusFilters');

disp('**********      buildNeuronsPostSpikeFilters      **********');
postSpikeHistory = buildNeuronsPostSpikeFilters(dt, lastPeakHistory, bForHistory, numOfBaseVectorsHistory, refreactoryPeriodArr, ISIPeakArr);
save('postSpikeHistory','postSpikeHistory');
% build coupling BaseVectors

disp('**********      buildBaseVectorsForPostSpikeAndCoupling      **********');
[~, ~, couplingBaseVectors] = buildBaseVectorsForPostSpikeAndCoupling(numOfBaseVectorsCoupling ,dt, [dt lastPeakCoupling], bForCoupling, 0);
 save('couplingBaseVectors', 'couplingBaseVectors');

% 
% figure();
% plot(couplingBaseVectors);drawnow;
disp('**********      getSpikeHistoryDataForNeurons      **********');
spikeHistoryData = getSpikeHistoryDataForNeurons(spikes, numOfBaseVectorsHistory, postSpikeHistory);
save('spikeHistoryData', 'spikeHistoryData', '-v7.3');

disp('**********      getCouplingDataForNeurons      **********');
couplingData = getCouplingDataForNeurons(spikes, numOfBaseVectorsCoupling, couplingBaseVectors);
save('couplingData', 'couplingData', '-v7.3');

end