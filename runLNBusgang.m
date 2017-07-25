function [expFunction, xData,yData]  = runLNBusgang(STA, nonRepeatStimulus, spikes, windowSize, numofParameters)

dataToRemove = 2000;
projectedStimulus = conv(nonRepeatStimulus, STA, 'same');
projectedStimulus = projectedStimulus(dataToRemove:end - dataToRemove);
spikes = spikes(dataToRemove:end - dataToRemove);
binnedprojectedStimulus = binData(projectedStimulus,windowSize);
binnedSpikes = binData(spikes,windowSize);
[expFunction, xData,yData] = buildNonLinFunctionFromData(binnedprojectedStimulus, binnedSpikes, numofParameters);
end