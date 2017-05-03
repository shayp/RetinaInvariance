datdir = './';  % directory where stimulus lives
load([datdir, 'Stim']);    % stimulus (temporal binary white noise)
load([datdir,'stimtimes']); % stim frame times in seconds (if desired)
load([datdir, 'SpTimes']); % load spike times (in units of stim frames)
ncells = length(SpTimes);  % number of neurons

Stim = Stim - mean(Stim);
% stimtimes(1:end) = stimtimes(1:end) / 10000;
% 
% % Pick a cell to work with
 cellnum = 1; % (1-2 are OFF cells; 3-4 are ON cells).
 tsp = SpTimes(cellnum).sp;
 wantedSampFactor = 20;
 [scaledSpikes, scaledStimulus] = changeSpikesAndStimulusRsolution(tsp, Stim,stimtimes, wantedSampFactor);
 lengthOfExp = length(scaledSpikes);
 trainfrac = .7;  % fraction of data to use for training
 ntrain = ceil(lengthOfExp*trainfrac);  % number of training samples
 ntest = lengthOfExp - ntrain; % number of test samples
 iitest = 1:ntest; % time indices for test
 iitrain = ntest+1:lengthOfExp;   % time indices for training
 stimtrain = scaledStimulus(iitrain,:); % training stimulus
 stimtest = scaledStimulus(iitest,:); % test stimulus
 spstrain = scaledSpikes(iitrain,:);
 spstest =  scaledSpikes(iitest,:);
 filterSizeBeforeSpike = 200;
 numOfBaseVectors = 10;
 [postSpiketimeVector,postSpikeBaseVectors, ~] = buildBaseVectorsForPostSpikeAndCoupling(numOfBaseVectors,0.05, [.1 2], .5);
 stimulusDesignMatrix = buildStimulusDesignMatrix(filterSizeBeforeSpike, stimtrain);
 filterSizeBeforeSpike = filterSizeBeforeSpike + 1;
 cellSTA = calculateSTA(stimulusDesignMatrix,spstrain);
 cellSTA(1)  = 0;

 spikeHistoryDesignMatrix = buildSpikeHistoryDesignMatrix(numOfBaseVectors, postSpikeBaseVectors, spstrain);

dataForLearnning.stimulusFiltetSize = filterSizeBeforeSpike;
dataForLearnning.binSizeInSecond = 1 / 500 ;
dataForLearnning.stimulusDesignMatrix = stimulusDesignMatrix;
dataForLearnning.spikeHistoryDesignMatrix = spikeHistoryDesignMatrix';
dataForLearnning.spikesccurence = find(spstrain);
dataForLearnning.dataLen = length(spstrain);
dataForLearnning.spikesTrain = spstrain;
% 
  opts = optimset('Gradobj','on','Hessian','off','display','iter-detailed');

 learnedParameters = [cellSTA' 0 0 0 0 0 0 0 0 0 0];
 %learnedParameters = [cellSTA'];
  lossfun = @(learnedParameters)Loss_GLM_logli_exp(learnedParameters,dataForLearnning); % set negative log-likelihood as loss func
  filtML = fminunc(lossfun,learnedParameters,opts);
