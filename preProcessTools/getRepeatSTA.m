function repSTA = getRepeatSTA(Stimulus,dataRepeats, stimulusFilterLength)
numOfSikes = sum(sum(dataRepeats));
repSTA = zeros(stimulusFilterLength,1);
numOfRepeats = size(dataRepeats, 1);
for i  = 1:numOfRepeats
    repSTA = repSTA + (dataRepeats(i,:) * Stimulus)';
end
repSTA = repSTA / numOfSikes;

end