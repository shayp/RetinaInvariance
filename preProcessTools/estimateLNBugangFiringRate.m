function lnBusgangFiringRate = estimateLNBugangFiringRate(stimulus, STA, expFunction, windowSize, deltaT)

projectedStimulus = conv(stimulus, STA, 'same');
binnedprojectedStimulus = binData(projectedStimulus,windowSize);
lnBusgangFiringRate = expFunction(binnedprojectedStimulus) * deltaT;
end