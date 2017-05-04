%%
datdir = './';  % directory where stimulus lives
load([datdir, 'Stim']);    % stimulus (temporal binary white noise)
load([datdir,'stimtimes']); % stim frame times in seconds (if desired)
load([datdir, 'SpTimes']); % load spike times (in units of stim frames)
ncells = length(SpTimes);  % number of neurons


% % Pick a cell to work with
 cellnum = 1; % (1-2 are OFF cells; 3-4 are ON cells).
 
 tsp = SpTimes(cellnum).sp;
 
 % Define the wanted factor to work with dseired tempral resplution(2ms)
 wantedSampFactor = 20;
 
 % We change the resolutin of the spikes and stimulus
 [scaledSpikes, scaledStimulus] = changeSpikesAndStimulusRsolution(tsp, Stim,stimtimes, wantedSampFactor);
 
 % We spilt the data for train and test
 lengthOfExp = length(scaledSpikes);
 trainfrac = .3;  % fraction of data to use for training
 ntrain = ceil(lengthOfExp*trainfrac);  % number of training samples
 ntest = lengthOfExp - ntrain; % number of test samples
 iitest = 1:ntest; % time indices for test
 iitrain = ntest+1:lengthOfExp;   % time indices for training
 stimtrain = scaledStimulus(iitrain,:); % training stimulus
 stimtest = scaledStimulus(iitest,:); % test stimulus
 spstrain = scaledSpikes(iitrain,:);
 spstest =  scaledSpikes(iitest,:);
 
 fprintf('Taining: %d spikes', sum(spstrain));
 fprintf('Testing: %d spikes', sum(spstest));
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
 teststimulusDesignMatrix = buildStimulusDesignMatrix(filterSizeBeforeSpike, stimtest);
 cellSTA = calculateSTA(stimulusDesignMatrix,spstrain);

 spikeHistoryDesignMatrix = buildSpikeHistoryDesignMatrix(numOfBaseVectors, postSpikeBaseVectors, spstrain);
 testspikeHistoryDesignMatrix = buildSpikeHistoryDesignMatrix(numOfBaseVectors, postSpikeBaseVectors, spstest);
 
dataForLearnning.stimulusFiltetSize = filterSizeBeforeSpike;
dataForLearnning.binSizeInSecond = 1 / 500 ;
dataForLearnning.stimulusDesignMatrix = stimulusDesignMatrix;
dataForLearnning.spikeHistoryDesignMatrix = spikeHistoryDesignMatrix';
dataForLearnning.spikesccurence = find(spstrain);
dataForLearnning.dataLen = length(spstrain);
dataForLearnning.spikesTrain = spstrain;

dataForTesting.stimulusFiltetSize = filterSizeBeforeSpike;
dataForTesting.binSizeInSecond = 1 / 500 ;
dataForTesting.stimulusDesignMatrix = teststimulusDesignMatrix;
dataForTesting.spikeHistoryDesignMatrix = testspikeHistoryDesignMatrix';
dataForTesting.spikesccurence = find(spstest);
dataForTesting.dataLen = length(spstest);
dataForTesting.spikesTrain = spstest;
% 
  opts = optimset('Gradobj','on','Hessian','on','display','iter-detailed');
  cellPostSpike = zeros(1,numOfBaseVectors);
  learnedParameters = [cellSTA' cellPostSpike];
  %learnedParameters = cellSTA';

% This matrix computes differences between adjacent coeffs
Dx1 = spdiags(ones(filterSizeBeforeSpike - 1,1)*[-1 1],0:1,filterSizeBeforeSpike-2,filterSizeBeforeSpike - 1);
Dx = Dx1'*Dx1; % computes squared diffs
% Select smoothing penalty by cross-validation 
lamvals = 2.^(1:20); % grid of lambda values (ridge parameters)
nlam = length(lamvals);

% Embed Dx matrix in matrix with one extra row/column for constant coeff
D = blkdiag(0,Dx); 

% Allocate space for train and test errors
negLtrain_sm = zeros(nlam,1);  % training error
negLogTest = zeros(nlam,1);  % training error

%%
%Initialize first parametets
negLikelihhod = @(prs)Loss_LN_logli_exp(prs,dataForLearnning); % set negative log-likelihood as loss func
learnedParameters = fminunc(negLikelihhod,learnedParameters(1:filterSizeBeforeSpike),opts);
learnedParameters = [learnedParameters cellPostSpike]; 
w_smooth = zeros(length(learnedParameters),nlam); % filters for each lambda


%%
figure();
hold on;
plot(learnedParameters(1:filterSizeBeforeSpike));drawnow;

for jj = 1:nlam
    wmap = learnedParameters; 
    % Compute MAP estimate
    Cinv = lamvals(jj)*D; % set inverse prior covariance
    negLikelihhod = @(prs)Loss_GLM_logli_exp(prs,dataForLearnning); % set negative log-likelihood as loss func
    lossfun = @(prs)glmneglogposterior(prs,negLikelihhod,Cinv, filterSizeBeforeSpike);
    learnedParameters = fminunc(lossfun,wmap,opts);
    w_smooth(:,jj) = learnedParameters;
    negLtrain_sm(jj) = Loss_GLM_logli_exp(learnedParameters, dataForLearnning);
    negLogTest(jj) = Loss_GLM_logli_exp(learnedParameters, dataForTesting);
    plot(learnedParameters(1:filterSizeBeforeSpike));
    xlabel('time before spike');drawnow; pause(.5);
end
hold off;
  figure();

  [~,imin] = min(negLogTest);
    subplot(4,1,1);
    plot(w_smooth(1:filterSizeBeforeSpike,imin));
    hold on;
    plot(cellSTA(1:end) - mean(cellSTA));
  spikeHistoryVector = w_smooth(end - numOfBaseVectors + 1 :end,imin)' * postSpikeBaseVectors';
  subplot(4,1,2);
  plot(spikeHistoryVector);
  subplot(4,1,3);
  plot(-negLtrain_sm);
  subplot(4,1,4);
  plot(-negLogTest);
