function [result_GLM_Full, result_GLM_Partial, result_LN,deltaT, cellSTA] = runGLM(neuronIndex, Stim, stimtimes, SpTimes, couplenNeurons)
%% Initlaization
addpath('./LossFunctions', './preProcessTools');
% set spikes var
tsp = SpTimes(neuronIndex).sp;
binsInSecond = 500;
deltaT = 1 / binsInSecond;
filterSizeBeforeSpike = 200;

numOfCoupledNeurons = length(couplenNeurons);
minLastSpike = tsp(length(tsp));

for i = 1:numOfCoupledNeurons
    % We p the coupled neurons in other struct for future use
    neuron(i).cellNumber = couplenNeurons(i); 
    neuron(i).spikeTimes = SpTimes(couplenNeurons(i)).sp;
    
    % We want to build a binary array of spikes, in the size of the minmal last
    % spike in all neurons that we use in current GLM Run
    if neuron(i).spikeTimes(length(neuron(i).spikeTimes)) < minLastSpike
        minLastSpike = neuron(i).spikeTimes(length(neuron(i).spikeTimes));
    end
end

% Define the wanted factor to work with dseired tempral resplution(2ms)
wantedSampFactor = 20;

% We change the resolutin of the spikes and stimulus
[scaledSpikes, scaledStimulus, rawSpikesVector,stimulusSampleVector, lastStimulus] = changeSpikesAndStimulusRsolution(tsp, Stim, stimtimes, wantedSampFactor, minLastSpike);


[neuronRawSpikes, neuronsSclaedSpikes] = getNeuronsRawSpikesSeries(numOfCoupledNeurons, neuron, lastStimulus, wantedSampFactor, length(scaledSpikes));
lengthOfExpRaw = length(rawSpikesVector);

% We spilt the data for train and test
lengthOfExp = length(scaledSpikes);
 
% fraction of data to use for training
trainfrac = .8;  
 
% number of training samples
ntrain = ceil(lengthOfExp*trainfrac);  

% number of test samples
ntest = lengthOfExp - ntrain; 

% time indices for training
iitrain = 1:ntrain;  

% time indices for test
iitest = ntrain + 1:lengthOfExp;
 
% training stimulus
stimtrain = scaledStimulus(iitrain);
% test stimulus
stimtest = scaledStimulus(iitest);

% Train spikes
spstrain = scaledSpikes(iitrain);
% Test spikes
spstest =  scaledSpikes(iitest);
coupledTrain = neuronsSclaedSpikes(:,iitrain);
coupledTest = neuronsSclaedSpikes(:,iitest);

% Print num of spikes 
fprintf('Taining scaled: %d spikes\n', sum(spstrain));
fprintf('Testing scaled: %d spikes\n', sum(spstest));

 
%% Post spike base vectors

% Define number of base vectors for post spike filter
numOfBaseVectors = 4;
lastPeak = 0.05;
dt = 0.001;
hpeaks = [0.001 lastPeak];
b = 0.005;
 
% build post spike BaseVectors
[postSpiketimeVector,postSpikeBaseVectors, originalBaseVectors] = buildBaseVectorsForPostSpikeAndCoupling(numOfBaseVectors,dt,hpeaks, b);
% Update the size after base vectors build(Can be changed)
numOfBaseVectors = size(postSpikeBaseVectors,2);
resizeBaseVectors = imresize(postSpikeBaseVectors, [40 numOfBaseVectors]);
postSpikeBaseVectors = resizeBaseVectors;
%% Design Matrix build

% We build the stimulus design matrix for train data
trainStimulusDesignMatrix = buildStimulusDesignMatrix(filterSizeBeforeSpike, stimtrain);

% We build the stimulus design matrix for test data
testStimulusDesignMatrix = buildStimulusDesignMatrix(filterSizeBeforeSpike, stimtest);

% We calculate the cell STA(With no zero-mean)
cellSTA = calculateSTA(trainStimulusDesignMatrix,spstrain);
cellSTA(1) = 0;

% We build spike history design matrix
[trainSpikeHistoryDesignMatrix, testSpikeHistoryDesignMatrix] = buildSpikeHistoryDesignMatrix(numOfBaseVectors, numOfCoupledNeurons,...
    postSpikeBaseVectors, length(spstrain),length(spstest), coupledTrain, coupledTest);

%% Optimization problem params init

% We create struct in order to gather data for the optimization
% problem(tesing and training)
dataForLearnning.stimulusFiltetSize = filterSizeBeforeSpike;
dataForLearnning.binSizeInSecond = 1 / binsInSecond ;
dataForLearnning.stimulusDesignMatrix = trainStimulusDesignMatrix;
dataForLearnning.spikeHistoryDesignMatrix = trainSpikeHistoryDesignMatrix';
dataForLearnning.dataLen = length(spstrain);
dataForLearnning.spikesTrain = spstrain;

dataForTesting.stimulusFiltetSize = filterSizeBeforeSpike;
dataForTesting.binSizeInSecond = 1 / binsInSecond ;
dataForTesting.stimulusDesignMatrix = testStimulusDesignMatrix;
dataForTesting.spikeHistoryDesignMatrix = testSpikeHistoryDesignMatrix';
dataForTesting.dataLen = length(spstest);
dataForTesting.spikesTrain = spstest;

% set options for optimizatiom problem
opts = optimset('Gradobj','on','Hessian','on','display','iter-detailed');

% Set start parameters for the optimization problem
cellPostSpike =  zeros(1,numOfBaseVectors * numOfCoupledNeurons);
meanFiringRate = rand();

% This matrix computes differences between adjacent coeffs
Dx1 = spdiags(ones(filterSizeBeforeSpike,1)*[-1 1],0:1,filterSizeBeforeSpike - 1,filterSizeBeforeSpike - 1);

% computes squared diffs(in order to get (w1 - w2).^2 penalty
Dx = Dx1'*Dx1; 

% Select lambda smoothing penalty by cross-validation 
 % grid of lambda values (ridge parameters)
lambdavals = (2).^(6:17);
nlambda = length(lambdavals);

% Embed Dx matrix in matrix with one extra row/column for constant coeff
D = blkdiag(0,Dx); 

%% Run optimization problem with diffrent lambdas
params_LN_Bias = [cellSTA' meanFiringRate]';
params_GLM_Full = [cellSTA' cellPostSpike meanFiringRate]';
params_GLM_Partial = [cellPostSpike meanFiringRate];
params_LN_No_Bias = cellSTA;

learned_LN_Bias = zeros(length(params_LN_Bias),nlambda);
learned_GLM_Full = zeros(length(params_GLM_Full),nlambda);
learned_GLM_Partial = zeros(length(params_GLM_Full),nlambda);
learned_LN_No_Bias = zeros(length(params_LN_No_Bias),nlambda);

% Allocate space for train and test errors
LL_GLM_Full_Train = zeros(nlambda,1);
LL_GLM_Full_Test = zeros(nlambda,1); 
LL_GLM_Partial_Train = zeros(nlambda,1);
LL_GLM_Partial_Test = zeros(nlambda,1);
LL_LN_Train = zeros(nlambda,1);
LL_LN_Test = zeros(nlambda,1);

fig = figure('visible', 'off');
hold on;
% Run for each lambda, and learn the the parameters
for i = 1:nlambda
    i
    % Compute the inverce covariance matrix with the current lambda
    Cinv = lambdavals(i) * D; % set inverse prior covariance
    
    % Learn  LN model with bias
    negLikelihhod = @(prs)Loss_LN_Bias(prs,dataForLearnning);
    lossfun = @(prs)Loss_posterior_LN_bias(prs,negLikelihhod,Cinv, filterSizeBeforeSpike);
    learned_LN_Bias(:,i) = fminunc(lossfun,params_LN_Bias,opts);
    LL_LN_Train(i) = Loss_LN_Bias(learned_LN_Bias(:,i), dataForLearnning);
    LL_LN_Test(i) = Loss_LN_Bias(learned_LN_Bias(:,i), dataForTesting);

    % Learn Full GLM
    negLikelihhod = @(prs)Loss_GLM_Full(prs,dataForLearnning);
    lossfun = @(prs)Loss_posterior_GLM(prs,negLikelihhod,Cinv, filterSizeBeforeSpike);
    learned_GLM_Full(:,i) = fminunc(lossfun,params_GLM_Full,opts);
    params_GLM_Full = learned_GLM_Full(:,i);
    LL_GLM_Full_Train(i) = Loss_GLM_Full(learned_GLM_Full(:,i), dataForLearnning);
    LL_GLM_Full_Test(i) = Loss_GLM_Full(learned_GLM_Full(:,i), dataForTesting);
    
    % Learn LN no bias
    negLikelihhod = @(prs)Loss_LN__NoBias(prs,dataForLearnning);
    lossfun = @(prs)Loss_posterior_LN_NoBias(prs,negLikelihhod,Cinv);
    learned_LN_No_Bias(:,i) = fminunc(lossfun,params_LN_No_Bias,opts);
    params_LN_No_Bias = learned_LN_No_Bias(:,i);
    
    % Learn GLM partial - using LN_No_Bias stimulus filter
    dataForLearnning.stimulusProjection = trainStimulusDesignMatrix * params_LN_No_Bias;
    dataForTesting.stimulusProjection = testStimulusDesignMatrix * params_LN_No_Bias;
    lossfun = @(prs)Loss_GLM_Partial_History(prs,dataForLearnning);
    params_GLM_Partial = fminunc(lossfun,params_GLM_Partial,opts);
    learned_GLM_Partial(:,i) = [params_LN_No_Bias' params_GLM_Partial]';
    LL_GLM_Partial_Train(i) = Loss_GLM_Partial_History(params_GLM_Partial, dataForLearnning);
    LL_GLM_Partial_Test(i) = Loss_GLM_Partial_History(params_GLM_Partial, dataForTesting);
end

% Get the minimum log likelihood index
[~,imin_LN_Bias] = min(LL_LN_Test);
[~,imin_GLM_Full] = min(LL_GLM_Full_Test);
[~,imin_GLM_Partial] = min(LL_GLM_Partial_Test);

bestParams_LN = learned_LN_Bias(:, imin_LN_Bias);
bestParams_GLM_Full = learned_GLM_Full(:, imin_GLM_Full);
bestParams_GLM_Partial = learned_GLM_Partial(:, imin_GLM_Partial);

for Index = 1:numOfCoupledNeurons
    result_GLM_Full.couplingFilters(Index,:) = bestParams_GLM_Full(filterSizeBeforeSpike + (Index  - 1) * numOfBaseVectors + 1 :filterSizeBeforeSpike + (Index) * numOfBaseVectors)' * postSpikeBaseVectors';
    result_GLM_Partial.couplingFilters(Index,:) = bestParams_GLM_Partial(filterSizeBeforeSpike + (Index  - 1) * numOfBaseVectors + 1 :filterSizeBeforeSpike + (Index) * numOfBaseVectors)' * postSpikeBaseVectors';
end
result_GLM_Full.StimulusFilter = bestParams_GLM_Full(1:filterSizeBeforeSpike);
result_GLM_Full.meanFiringRate = bestParams_GLM_Full(end);

result_GLM_Partial.StimulusFilter = bestParams_GLM_Partial(1:filterSizeBeforeSpike);
result_GLM_Partial.meanFiringRate = bestParams_GLM_Partial(end);

result_LN.StimulusFilter = bestParams_LN(1:filterSizeBeforeSpike);
result_LN.meanFiringRate = bestParams_LN(end);
%% Plot learned estimators
% 
% fig = figure('visible', 'off');
% 
% % Train likelihood
% subplot(3,1,1);
% plot(-negLogTrain);
% title('train likelihood');
% xlabel('lambda factor');
% ylabel('log likelihood');
% 
% % Test likelihood
% subplot(3,1,2);
% plot(-negLogTest);
% title('test likelihood');
% xlabel('lambda factor');
% ylabel('log likelihood');
% 
% % Test likelihood
% subplot(3,1,3);
% plot(meanFiringRateArray);
% title('Mean firing rate');
% xlabel('lambda factor');
% ylabel('mean value');
% savefig(fig,['./Graphs/Neuron_' num2str(neuronIndex) '_Likelihood']);

end
