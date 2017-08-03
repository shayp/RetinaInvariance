load SpikeTimeGaussian
load GaussianMovieData

SetGlobalVals();
load('Global');
% WW(:,1) - Event marker, marks important event on the screen
% StimTime - Array that contains the start time of a new event(can be
% reapeat or non reapet event)

% number of neurons to include
N = 251;
lengthOfSTAFilter = 5000;
numOfSessions = 6;
wantedSeries = 6;
numOfClusters = 6;
% define array for the sta
STAAll=zeros(lengthOfSTAFilter,N);

spikeQuality = ([TT(:).Quality]);
[~, sortedQualityInd] = sort(spikeQuality);
sortedQualityInd = fliplr(sortedQualityInd);
selectedenuron =144;
save (['sortedQualityInd'],'sortedQualityInd');  


% we find only the imortant indexes in the event marker
importantStimIndexes = find(WW(:,1)>150);
problematicIndexes = find(diff(importantStimIndexes) < 10) + 1;
importantStimIndexes(problematicIndexes) = [];
% Find first and last indexes for the wanted seesion
wantedSeriesFirstIndex = (globalVal.numberOfNonRepeatStimulus + globalVal.numberOfRepeatStimulus) * (wantedSeries - 1) + 1;
wantedSeriesLastIndex = (globalVal.numberOfNonRepeatStimulus + globalVal.numberOfRepeatStimulus) * (wantedSeries) + 1;


% Find start and end time for stimulus in the wanted session
nonRepeatStimulusStartTime = StimTime(wantedSeriesFirstIndex);
RepeatStimulusStartTime = StimTime(wantedSeriesFirstIndex + globalVal.numberOfNonRepeatStimulus);
endOfStimlus = StimTime(wantedSeriesLastIndex);
repeatStimulusTimes = StimTime((wantedSeriesFirstIndex + globalVal.numberOfNonRepeatStimulus: wantedSeriesLastIndex));
repeatStimulusTimes = repeatStimulusTimes - repeatStimulusTimes(1) + 1;

% Define the  repeat stimulus array
Stimulus = zeros(2,endOfStimlus - nonRepeatStimulusStartTime); 

%Remove mean from the stimulus(zero mean)
%WW(:,3) = WW(:,3) - mean(WW(:,3));

% Insert the non repeat stimulus
for j = wantedSeriesFirstIndex: wantedSeriesFirstIndex + globalVal.numberOfNonRepeatStimulus - 1
    % find the length for resize the stimuli in this session
    lengthToResize = StimTime(j + 1) - StimTime(j);
    % Get current stimulus from the movie data
    currentStimulus = WW(importantStimIndexes(j):importantStimIndexes(j + 1) - 1,3);
    currentStimulus = currentStimulus - mean(currentStimulus);
    % Extend the stimulus to the expiriment size
    extendStimulus = imresize(currentStimulus, [lengthToResize 1], 'nearest');
    
    % Define array for the meanning od each stimulus
    stimuliValArray = zeros(1,lengthToResize);
    
    stimuliValArray(1:lengthToResize) =  0;

    % Set the stimulus and the type of the stimlus
    Stimulus(1, StimTime(j) - StimTime(wantedSeriesFirstIndex) + 1:StimTime(j + 1) - StimTime(wantedSeriesFirstIndex)) = extendStimulus;
    Stimulus(2, StimTime(j) - StimTime(wantedSeriesFirstIndex) + 1:StimTime(j + 1) - StimTime(wantedSeriesFirstIndex)) = stimuliValArray;
end
disp('********************************');


% Run for each "small" session of stimulus
for j = wantedSeriesFirstIndex + globalVal.numberOfNonRepeatStimulus:wantedSeriesLastIndex -1
    
    % find the length for resize the stimuli in this session
    lengthToResize = StimTime(j + 1) - StimTime(j);
    
	currentStimulus = WW(importantStimIndexes(j) + 1:importantStimIndexes(j + 1),3);
    currentStimulus = currentStimulus - mean(currentStimulus);
    extendStimulus = zeros(1, lengthToResize);
    currentIndex = 1;
    for i = 1:length(currentStimulus)
        if mod(i, 3) == 0
            currentStep = 334;
        else
            currentStep = 333;
        end
        if currentIndex + currentStep - 1 >= lengthToResize
            break;
        end
        extendStimulus(currentIndex: currentIndex + currentStep - 1) = currentStimulus(i);
        currentIndex = currentIndex + currentStep;

    end
    
    extendStimulus(currentIndex:end) = currentStimulus(end);
    % Extend the stimulus to the expiriment size
%     extendStimulus = imresize(currentStimulus, [lengthToResize 1], 'nearest');
%     diff(find(diff(extendStimulus) ~= 0))
%     immse(compareStimulus, currentStimulus)
    % Define array for the meanning od each stimulus
    stimuliValArray = zeros(1,lengthToResize);
    stimuliValArray =  j - wantedSeriesFirstIndex - globalVal.numberOfNonRepeatStimulus + 1;
    
    % Set the stimulus and the type of the stimlus
    Stimulus(1, StimTime(j) - StimTime(wantedSeriesFirstIndex) + 1:StimTime(j + 1) - StimTime(wantedSeriesFirstIndex)) = extendStimulus;
    Stimulus(2, StimTime(j) - StimTime(wantedSeriesFirstIndex) + 1:StimTime(j + 1) - StimTime(wantedSeriesFirstIndex)) = stimuliValArray;
end

% Now we want to take the neurons and remove spikes before and after the
% expiriment time
 maxSpikesNum = 0;

  nonRepStimulusExtended = Stimulus(1,find(~Stimulus(2,:)));
  RepStimulusExtended = Stimulus(:,find(Stimulus(2,:)));
  nonRepstimulus = [nonRepStimulusExtended(1), nonRepStimulusExtended(find(diff(nonRepStimulusExtended(:)) ~= 0) + 1)];
  x = ones(1,1);
  nonRepstimulusTimes = [x, [find(diff(nonRepStimulusExtended(:)) ~= 0) + 1]'];
  Stim = nonRepstimulus;
  stimtimes = nonRepstimulusTimes;
  save ('NonRepeatStim','Stim');    
  save ('NonRepeatstimtimes','stimtimes');    
  save ('repeatStimulusTimes', 'repeatStimulusTimes');    
  save ('RepStimulusExtended', 'RepStimulusExtended');    
%%
NonRepTT = TT;
 for i = 1:N
    NonRepTT(i).sp = NonRepTT(i).sp(find(NonRepTT(i).sp > nonRepeatStimulusStartTime + lengthOfSTAFilter));
    NonRepTT(i).sp = NonRepTT(i).sp(find(NonRepTT(i).sp < RepeatStimulusStartTime));
    NonRepTT(i).sp = NonRepTT(i).sp - StimTime(wantedSeriesFirstIndex) + 1;
 end

SpTimes = NonRepTT;
  save (['SpTimes'],'SpTimes');    
  
  RepTT = TT;
  
for i = 1:N
    RepTT(i).sp = RepTT(i).sp(find(RepTT(i).sp > RepeatStimulusStartTime + lengthOfSTAFilter));
    RepTT(i).sp = RepTT(i).sp(find(RepTT(i).sp < endOfStimlus));
   
    RepTT(i).sp = RepTT(i).sp - StimTime(wantedSeriesFirstIndex + globalVal.numberOfNonRepeatStimulus) + 1;
 end

RepSpTimes = RepTT;
  save (['RepSpTimes'],'RepSpTimes');    
%%
RepSTA = zeros(1, lengthOfSTAFilter);
nonRepSTA = zeros(1, lengthOfSTAFilter);
length(RepTT(selectedenuron).sp)
length(NonRepTT(selectedenuron).sp)

RepTT(selectedenuron).sp = RepTT(selectedenuron).sp + StimTime(wantedSeriesFirstIndex + globalVal.numberOfNonRepeatStimulus) - StimTime(wantedSeriesFirstIndex);
repDesignMatrix = zeros(lengthOfSTAFilter,5000);
nonRepDesignMatrix = zeros(lengthOfSTAFilter,5000);
for i=selectedenuron
    for spikeIndex = 1:length(RepTT(i).sp)
        repDesignMatrix(1:lengthOfSTAFilter,spikeIndex) = Stimulus(1, RepTT(i).sp(spikeIndex) - lengthOfSTAFilter + 1:RepTT(i).sp(spikeIndex));
    end
    for spikeIndex = 1:length(NonRepTT(i).sp)
        nonRepDesignMatrix(1:lengthOfSTAFilter,spikeIndex) = Stimulus(1, NonRepTT(i).sp(spikeIndex) - lengthOfSTAFilter + 1:NonRepTT(i).sp(spikeIndex));
    end
    
   RepSTA(1,1:lengthOfSTAFilter) = sum(repDesignMatrix( :, :),2) / length(RepTT(i).sp);
   nonRepSTA(1,1:lengthOfSTAFilter) = sum(nonRepDesignMatrix(:, :),2) / length(NonRepTT(i).sp);
end

 figure();
 plot(1:lengthOfSTAFilter, RepSTA(1,:), 1:lengthOfSTAFilter,  nonRepSTA(1,:));
 legend('Rep', 'Non Rep');
