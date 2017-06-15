load SpikeTimeGaussian
load GaussianMovieData

SetGlobalVals();
load('Global');
% WW(:,1) - Event marker, marks important event on the screen
% StimTime - Array that contains the start time of a new event(can be
% reapeat or non reapet event)

% number of neurons to include
N = 251;
lengthOfSTAFilter = 4000;
wantedSeries = 6;
numOfClusters = 6;
selectedenuron = 33;
% define array for the sta
STAAll=zeros(lengthOfSTAFilter,N);

spikeQuality = ([TT(:).Quality]);
[~, sortedQualityInd] = sort(spikeQuality);
save (['sortedQualityInd'],'sortedQualityInd');    
% we find only the imortant indexes in the event marker
importantStimIndexes = find(WW(:,1)>150);
importantStimIndexes = importantStimIndexes(find(diff(importantStimIndexes)>10));
% Find first and last indexes for the wanted seesion
wantedSeriesFirstIndex = (globalVal.numberOfNonRepeatStimulus + globalVal.numberOfRepeatStimulus) * (wantedSeries - 1) + 1;
wantedSeriesLastIndex = (globalVal.numberOfNonRepeatStimulus + globalVal.numberOfRepeatStimulus) * (wantedSeries);

% Find start and end time for stimulus in the wanted session
nonRepeatStimulusStartTime = StimTime(wantedSeriesFirstIndex);
RepeatStimulusStartTime = StimTime(wantedSeriesFirstIndex + globalVal.numberOfNonRepeatStimulus);
endOfStimlus = StimTime(wantedSeriesLastIndex);
repeatStimulusTimes = StimTime((wantedSeriesFirstIndex + globalVal.numberOfNonRepeatStimulus: wantedSeriesLastIndex));
repeatStimulusTimes = repeatStimulusTimes - repeatStimulusTimes(1) + 1;
% define where the repeat stimulus starts
repeatStimulusStartInVector = RepeatStimulusStartTime - nonRepeatStimulusStartTime;

% Get only the important indexes of wanted session from the event marker
wantedseriesWWIndexes = importantStimIndexes(wantedSeriesFirstIndex:wantedSeriesLastIndex);

% Define the stimulus array
Stimulus = zeros(2,endOfStimlus - nonRepeatStimulusStartTime); 

%Remove mean from the stimulus(zero mean)
WW(:,3) = WW(:,3) - mean(WW(:,3));
for j = wantedSeriesFirstIndex: wantedSeriesFirstIndex + globalVal.numberOfNonRepeatStimulus - 1
    % find the length for resize the stimuli in this session
    lengthToResize = StimTime(j + 1) - StimTime(j);
    % Get current stimulus from the movie data
    currentStimulus = WW(importantStimIndexes(j):importantStimIndexes(j + 1) - 1,3);
    
     % Extend the stimulus to the expiriment size
    extendStimulus = imresize(currentStimulus, [lengthToResize 1], 'nearest');
    
    % Define array for the meanning od each stimulus
    stimuliValArray = zeros(1,lengthToResize);
    
    stimuliValArray(1:lengthToResize) =  0;

    
    % Set the stimulus and the type of the stimlus
    Stimulus(1, StimTime(j) - StimTime(wantedSeriesFirstIndex) + 1:StimTime(j + 1) - StimTime(wantedSeriesFirstIndex)) = extendStimulus;
    Stimulus(2, StimTime(j) - StimTime(wantedSeriesFirstIndex) + 1:StimTime(j + 1) - StimTime(wantedSeriesFirstIndex)) = stimuliValArray;
end
% Run for each "small" session of stimulus
for j = wantedSeriesFirstIndex + globalVal.numberOfNonRepeatStimulus:wantedSeriesLastIndex - 1
    
    % find the length for resize the stimuli in this session
    lengthToResize = StimTime(j + 1) - StimTime(j);
    
	currentStimulus = WW(importantStimIndexes(j) + 1:importantStimIndexes(j + 1) ,3);

    j
   lengthToResize
   currentStimulus(1)
   currentStimulus(end)

    % Extend the stimulus to the expiriment size
    extendStimulus = imresize(currentStimulus, [lengthToResize 1], 'nearest');
    
    % Define array for the meanning od each stimulus
    stimuliValArray = zeros(1,lengthToResize);
    
    stimuliValArray(1:lengthToResize) =  j - wantedSeriesFirstIndex - globalVal.numberOfNonRepeatStimulus + 1;
    
    % Set the stimulus and the type of the stimlus
    Stimulus(1, StimTime(j) - StimTime(wantedSeriesFirstIndex) + 1:StimTime(j + 1) - StimTime(wantedSeriesFirstIndex)) = extendStimulus;
    Stimulus(2, StimTime(j) - StimTime(wantedSeriesFirstIndex) + 1:StimTime(j + 1) - StimTime(wantedSeriesFirstIndex)) = stimuliValArray;
end

% Now we want to take the neurons and remove spikes before and after the
% expiriment time
 maxSpikesNum = 0;

  nonRepStimulusExtended = Stimulus(1,find(Stimulus(2,:) == 0));
  RepStimulusExtended = Stimulus(:,find(Stimulus(2,:) ~= 0));
  nonRepstimulus = [nonRepStimulusExtended(1), nonRepStimulusExtended(find(diff(nonRepStimulusExtended(:)) ~= 0) + 1)];
  x = ones(1,1);
  nonRepstimulusTimes = [x, [find(diff(nonRepStimulusExtended(:)) ~= 0) + 1]'];
  Stim = nonRepstimulus;
  stimtimes = nonRepstimulusTimes;
  save (['Stim'],'Stim');    
  save (['stimtimes'],'stimtimes');    
  save (['repeatStimulusTimes'],'repeatStimulusTimes');    
  save (['RepStimulusExtended'],'RepStimulusExtended');    
%%
NonRepTT = TT;
 for i = 1:N
    NonRepTT(i).sp = NonRepTT(i).sp(find(NonRepTT(i).sp > nonRepeatStimulusStartTime + lengthOfSTAFilter));
    NonRepTT(i).sp = NonRepTT(i).sp(find(NonRepTT(i).sp < RepeatStimulusStartTime));
    if length(NonRepTT(i).sp) > maxSpikesNum
        maxSpikesNum = length(NonRepTT(i).sp);
    end
    
    NonRepTT(i).sp = NonRepTT(i).sp - StimTime(wantedSeriesFirstIndex) + 1;
 end

SpTimes = NonRepTT;
  save (['SpTimes'],'SpTimes');    
  
  RepTT = TT;
  
   for i = 1:N
    RepTT(i).sp = RepTT(i).sp(find(RepTT(i).sp > RepeatStimulusStartTime + lengthOfSTAFilter));
    RepTT(i).sp = RepTT(i).sp(find(RepTT(i).sp < endOfStimlus));
    if length(RepTT(i).sp) > maxSpikesNum
        maxSpikesNum = length(RepTT(i).sp);
    end
    
    RepTT(i).sp = RepTT(i).sp - StimTime(wantedSeriesFirstIndex + globalVal.numberOfNonRepeatStimulus) + 1;
 end

RepSpTimes = RepTT;
  save (['RepSpTimes'],'RepSpTimes');    
%%
RepSTA = zeros(N, lengthOfSTAFilter);
nonRepSTA = zeros(N, lengthOfSTAFilter);

RepTT(selectedenuron).sp = RepTT(selectedenuron).sp + StimTime(wantedSeriesFirstIndex + globalVal.numberOfNonRepeatStimulus) - StimTime(wantedSeriesFirstIndex);
repDesignMatrix = zeros(1,lengthOfSTAFilter,5000);
nonRepDesignMatrix = zeros(1,lengthOfSTAFilter,5000);
for i=selectedenuron
    for spikeIndex = 1:length(RepTT(i).sp);
        repDesignMatrix(i,1:lengthOfSTAFilter,spikeIndex) = Stimulus(1, RepTT(i).sp(spikeIndex) - lengthOfSTAFilter:RepTT(i).sp(spikeIndex) - 1);
    end
    for spikeIndex = 1:length(NonRepTT(i).sp);
        nonRepDesignMatrix(i,1:lengthOfSTAFilter,spikeIndex) = Stimulus(1, NonRepTT(i).sp(spikeIndex) - lengthOfSTAFilter:NonRepTT(i).sp(spikeIndex) - 1);
    end
    
   RepSTA(i,1:lengthOfSTAFilter) = sum(repDesignMatrix(i, :, :),3) / length(RepTT(i).sp);
   nonRepSTA(i,1:lengthOfSTAFilter) = sum(nonRepDesignMatrix(i, :, :),3) / length(NonRepTT(i).sp);
end

%  for i=selectedenuron
%     for spikeIndex = 1:length(RepTT(i).sp);
%         RepSTA(i,:) = RepSTA(i,:) + Stimulus(1, RepTT(i).sp(spikeIndex) - lengthOfSTAFilter:RepTT(i).sp(spikeIndex) - 1);
%     end
%     RepSTA(i,:) = RepSTA(i,:) / length(RepTT(i).sp);
%    
%     for spikeIndex = 1:length(NonRepTT(i).sp);
%         nonRepSTA(i,:) = nonRepSTA(i,:) + Stimulus(1, NonRepTT(i).sp(spikeIndex) - lengthOfSTAFilter:NonRepTT(i).sp(spikeIndex) - 1);
%     end
%     length(NonRepTT(i).sp)
%     nonRepSTA(i,:) = nonRepSTA(i,:) / length(NonRepTT(i).sp);
%  end
 figure();
 plot(1:lengthOfSTAFilter, RepSTA(selectedenuron,:), 1:lengthOfSTAFilter,  nonRepSTA(selectedenuron,:));
 
% Removelist = [192,161,141,59];
% STA(Removelist(1:end),:) = [];
% T = clusterdata(STA, 'linkage', 'average','maxclust', numOfClusters);
% clusturedSTA = zeros(numOfClusters, lengthOfSTAFilter);
% clusterCount = zeros(numOfClusters,1);
% 
% ttk = (-lengthOfSTAFilter+1:0) / 10000;
% figure();
%  curveColors = char({'yellow', 'red', 'magenta','green','cyan', 'black'});
%  size(STA,1)
% for i=1:size(STA,1);
%     [i i + length(find(Removelist(:) < i))  T(i)]
%     hold on;
%     tmp = STA(i,:);
%     clusturedSTA(T(i),1:lengthOfSTAFilter) = clusturedSTA(T(i),1:lengthOfSTAFilter) + tmp;
%     clusterCount(T(i)) = clusterCount(T(i)) + 1;
%     subplot(2,1,1);
%     
%     title('Clustered neurons');
%     xlabel('Time Before Spilke (s)');
%     ylabel('STA');
%     plot(ttk,tmp,curveColors(T(i)),'LineWidth',1);
% 
%     hold off;
% end
% clusterCount
% for i=1:numOfClusters
%     hold on;
%     if clusterCount(i) ~= 0
%     clusturedSTA(i,:) = clusturedSTA(i,:) / clusterCount(i);
%     subplot(2,1,2);
%     title('Mean of Clustered neurons');
%     xlabel('Time Before Spilke (s)');
%     ylabel('STA');
%     plot(ttk,clusturedSTA(i,:),curveColors(i),'DisplayName',['cluster ' num2str(i) ' - ' num2str(clusterCount(i)) ' neurons'],'LineWidth',1);
%     end
%     hold off;
% end
% legend('show')
% 
%     postSpikeFilter = 150;
%     postSpikeAVg = zeros(N,postSpikeFilter + 1);
% for i=1:1
%     numOfSpikes = length(TT(i).sp);
%     currentSpikesVector = ismember( 1:TT(i).sp(numOfSpikes),TT(i).sp);
%     for j=1:numOfSpikes
%         if TT(i).sp(j) < TT(i).sp(numOfSpikes) - postSpikeFilter
% 
%             length(TT(i).sp(j) + 1:TT(i).sp(j)  + postSpikeFilter)
%             length(postSpikeAVg(i,2:end))
%             postSpikeAVg(i,2:end) = postSpikeAVg(i,2:end) + currentSpikesVector(TT(i).sp(j) + 1:TT(i).sp(j)  + postSpikeFilter);
%             postSpikeAVg(i,1) = postSpikeAVg(i,1) + 1;
%         end
%     end
%     %postSpikeAVg(i,2:end) = postSpikeAVg(i,2:end) / postSpikeAVg(i,1);
% end
% figure();
% plot(postSpikeAVg(1,2:end) / max(postSpikeAVg(1,2:end)));
