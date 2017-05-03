datdir = './';  % directory where stimulus lives
load([datdir, 'Stim']);    % stimulus (temporal binary white noise)
load([datdir,'stimtimes']); % stim frame times in seconds (if desired)
load([datdir, 'SpTimes']); % load spike times (in units of stim frames)
ncells = length(SpTimes);  % number of neurons

% Change stimulus to zero mean intensity
Stim = Stim - mean(Stim);

% % Pick a cell to work with
 cellnum = 1; % (1-2 are OFF cells; 3-4 are ON cells).
 
 tsp = SpTimes(cellnum).sp;
 
 % Define the wanted factor to work with dseired tempral resplution(2ms)
 wantedSampFactor = 20;
 
 % We change the resolutin of the spikes and stimulus
 [scaledSpikes, scaledStimulus] = changeSpikesAndStimulusRsolution(tsp, Stim,stimtimes, wantedSampFactor);
 
 % We spilt the data for train and test
 lengthOfExp = length(scaledSpikes);
 trainfrac = .05;  % fraction of data to use for training
 ntrain = ceil(lengthOfExp*trainfrac);  % number of training samples
 ntest = lengthOfExp - ntrain; % number of test samples
 iitest = 1:ntest; % time indices for test
 iitrain = ntest+1:lengthOfExp;   % time indices for training
 stimtrain = scaledStimulus(iitrain,:); % training stimulus
 stimtest = scaledStimulus(iitest,:); % test stimulus
 spstrain = scaledSpikes(iitrain,:);
 spstest =  scaledSpikes(iitest,:);
 
 % Define the wanted size of stimulus filter before spike
 filterSizeBeforeSpike = 200;
 
 numOfBaseVectors = 10;
 
 % Build base vectors for post spike history
 lastPeak = 4;
 dt = 1;
 [postSpiketimeVector,postSpikeBaseVectors, originalBaseVectors] = buildBaseVectorsForPostSpikeAndCoupling(numOfBaseVectors,dt, [0.05 lastPeak], lastPeak / 2);
 numOfBaseVectors = size(postSpikeBaseVectors,2);

 % We build the stimulus design matrix
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
  opts = optimset('Gradobj','on','Hessian','on','display','iter-detailed');
  cellPostSpike = zeros(1,numOfBaseVectors);
  learnedParameters = [cellSTA'];

  % Precompute some quantities (X'X and X'*y) for training and test data
  Imat = eye(filterSizeBeforeSpike); % identity matrix of size of filter + const
  Imat(1,1) = 0; % remove penalty on constant dc offset
  Cinv = 1024*Imat; % set inverse prior covariance
  size(Cinv)
  negLikelihhod = @(learnedParameters)Loss_GLM_logli_exp(learnedParameters,dataForLearnning); % set negative log-likelihood as loss func
  lossfun = @(prs)neglogposterior(learnedParameters,negLikelihhod,Cinv);
  filtML = fminunc(lossfun,learnedParameters,opts);
  plot(filtML);
  hold on;
  plot(cellSTA , '.');
