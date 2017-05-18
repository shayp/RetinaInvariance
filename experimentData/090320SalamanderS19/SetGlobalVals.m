function SetGlobalVals()
globalVal.numberOfLongExpiriments = 6;
globalVal.signBeforeSmallStimuliSeries = 180;
globalVal.signBeforeRepeatStimulus = 120;
globalVal.signForBigStimulusSeries = 255;
globalVal.lengthOfSTAFilter = 5000;
globalVal.numberOfNonRepeatStimulus = 11;
globalVal.numberOfRepeatStimulus = 10;
globalVal.lengthOfNonRepeatStimuli = 1800;
globalVal.lengthOfRepeatStimuli = 601;
%globalVal.lengthOfStimuliSession = numberOfNonRepeatStimulus * lengthOfNonRepeatStimuli + globalVal.numberOfRepeatStimuli * lengthOfRepeatStimuli;
globalVal.sistimulusValIndex = 1;
globalVal.SampleSize = 10000;
globalVal.screenRate = 30;
globalVal.stimulusSizeInSampleSize = globalVal.SampleSize / globalVal.screenRate;
save ('Global', 'globalVal');
end