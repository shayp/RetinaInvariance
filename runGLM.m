function [scaledStimulus, couplingFilters, learnedSTA, deltaT, meanFiringRate] = runGLM(neuronIndex, Stim, stimtimes, SpTimes, couplenNeurons)
%% Initlaization

% set spikes var
tsp = SpTimes(neuronIndex).sp;
%tsp = tsp(1:1500);
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
% bla(scaledStimulus,find(scaledSpikes), 200);
% bla(stimulusSampleVector,find(rawSpikesVector), 4000);

[neuronRawSpikes, neuronsSclaedSpikes] = getNeuronsRawSpikesSeries(numOfCoupledNeurons, neuron, lastStimulus, wantedSampFactor, length(scaledSpikes));
lengthOfExpRaw = length(rawSpikesVector);

% We spilt the data for train and test
lengthOfExp = length(scaledSpikes);
 
% fraction of data to use for training
trainfrac = .25;  
 
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
 
lastPeak = 0.060;
dt = 0.002;
hpeaks = [0.001 lastPeak];
b = 0.005;
 
% build post spike BaseVectors
[postSpiketimeVector,postSpikeBaseVectors, originalBaseVectors] = buildBaseVectorsForPostSpikeAndCoupling(numOfBaseVectors,dt,hpeaks, b);
% Update the size after base vectors build(Can be changed)
numOfBaseVectors = size(postSpikeBaseVectors,2);

% % % % % Plot base vectors
% figure();
% plot(postSpikeBaseVectors);
% title('Base vectors for post spike history');
% xlabel('Time after spike');drawnow;
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
cellPostSpike =  0.01 * rand(1,numOfBaseVectors * numOfCoupledNeurons);
meanFiringRate = 4;
learnedParameters = [cellSTA' cellPostSpike meanFiringRate];

% This matrix computes differences between adjacent coeffs
Dx1 = spdiags(ones(filterSizeBeforeSpike,1)*[-1 1],0:1,filterSizeBeforeSpike - 1,filterSizeBeforeSpike - 1);

% computes squared diffs(in order to get (w1 - w2).^2 penalty
Dx = Dx1'*Dx1; 

% Select lambda smoothing penalty by cross-validation 
 % grid of lambda values (ridge parameters)
lambdavals = (2).^(1:20);
nlambda = length(lambdavals);

% Embed Dx matrix in matrix with one extra row/column for constant coeff
D = blkdiag(0,Dx); 

% Allocate space for train and test errors
negLogTrain = zeros(nlambda,1); 
negLogTest = zeros(nlambda,1);  
lambdaLearrnedParameters = zeros(length(learnedParameters),nlambda);
meanFiringRateArray = zeros(1,nlambda);

%% Run optimization problem with diffrent lambdas

startParams = [cellSTA' cellPostSpike meanFiringRate];
learnedParameters = startParams;
figure();
hold on;

% Run for each lambda, and learn the the parameters
for i = 1:nlambda
    i
    % We plot the learned STA of current iteration
    plot(learnedParameters(1:filterSizeBeforeSpike));
    ylabel('intensity');
    title([num2str(i)]);
    xlabel('time before spike');drawnow; pause(.5);

    % The cirrent entry parameters are the previous learing estimatror
    currentEntryParameters = learnedParameters; 
    
    % Compute the inverce covariance matrix with the current lambda
    Cinv = lambdavals(i) * D; % set inverse prior covariance
    
    % We use 2 loss function, the first is negative log likelihood based on
    % inhomginious poission process. The second is smoothing the STA
    % estimatror using penalty on the difference between two linked weights.
    
    % The negative log likelihood function
    negLikelihhod = @(prs)Loss_GLM_logli_exp(prs,dataForLearnning);
    
    % The negative log posteriot funnction, for smoothing STA
    lossfun = @(prs)glmneglogposterior(prs,negLikelihhod,Cinv, filterSizeBeforeSpike);
    
    % Call The optimization problem solver
    learnedParameters = fminunc(lossfun,currentEntryParameters,opts);

    % We save the learned parameters
    lambdaLearrnedParameters(:,i) = learnedParameters;
    
    % We calculate the negative log likelihood value for the current
    % estimator
    negLogTrain(i) = Loss_GLM_logli_exp(learnedParameters, dataForLearnning);
    negLogTest(i) = Loss_GLM_logli_exp(learnedParameters, dataForTesting);
    meanFiringRateArray(i) = learnedParameters(end);
end

hold off;

% Get the minimum log likelihood index
[~,imin] = min(negLogTest);
choosedParams = lambdaLearrnedParameters(:,imin);

% Calculate the spike history vector based on the parameters learned
spikeHistoryVector = choosedParams(filterSizeBeforeSpike + 1 :filterSizeBeforeSpike + numOfBaseVectors)' * postSpikeBaseVectors';
couplingFilters = zeros(numOfCoupledNeurons, size(postSpikeBaseVectors,1));

for neuronIndex = 1:numOfCoupledNeurons
    couplingFilters(neuronIndex,:) = choosedParams(filterSizeBeforeSpike + (neuronIndex  - 1) * numOfBaseVectors + 1 :filterSizeBeforeSpike + (neuronIndex) * numOfBaseVectors)' * postSpikeBaseVectors';
end
learnedSTA = lambdaLearrnedParameters(1:filterSizeBeforeSpike,imin);
meanFiringRate = lambdaLearrnedParameters(end,imin);
%% Plot learned estimators

figure();

% Plot STA estimator
subplot(3,2,1);
plot(lambdaLearrnedParameters(1:filterSizeBeforeSpike,imin));
hold on;
plot(cellSTA(1:end));
legend('learned STA','Expiriment STA');
xlabel('Time before spike');
ylabel('intensity');
title('STA estimatror');

subplot(3,2,2);
plot(meanFiringRateArray);
title('mean firing rate');
xlabel('lambda factor');
ylabel('mean firing rate');

% Plot leaned spike history filter
subplot(3,2,3);
plot(couplingFilters(1,:));
title('coupling filter');
xlabel('Time after  spike');
ylabel('Firing factor');

% % Plot leaned spike history filter
% subplot(3,2,4);
% plot(couplingFilters(2,:));
% title('coupling filter');
% xlabel('Time after spike');
% ylabel('Firing factor');

% Train likelihood
subplot(3,2,5);
plot(-negLogTrain);
title('train likelihood');
xlabel('lambda factor');
ylabel('log likelihood');

% Test likelihood
subplot(3,2,6);
plot(-negLogTest);
title('test likelihood');
xlabel('lambda factor');
ylabel('log likelihood');
end
