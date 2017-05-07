function runGLM(neuronIndex, Stim, stimtimes, SpTimes, couplenNeurons)
%% Initlaization

% set spikes var
tsp = SpTimes(neuronIndex).sp;
 
binsInSecond = 500;

numOfCoupledNeurons = length(couplenNeurons);

for i = 1:numOfCoupledNeurons
    neuron(i).cellNumber = couplenNeurons(i); 
    neuron(i).spikeTimes = SpTimes(couplenNeurons(i)).sp;
end

% Define the wanted factor to work with dseired tempral resplution(2ms)
wantedSampFactor = 20;
 
% We change the resolutin of the spikes and stimulus
[scaledSpikes, scaledStimulus, rawSpikesVector, lastStimulus] = changeSpikesAndStimulusRsolution(tsp, Stim, stimtimes, wantedSampFactor);

neuronRawStimulus = getNeuronsRawStimulusSeires(numOfCoupledNeurons, neuron, lastStimulus);

lengthOfExpRaw = length(rawSpikesVector);

% We spilt the data for train and test
lengthOfExp = length(scaledSpikes);
 
% fraction of data to use for training
trainfrac = .3;  
 
% number of training samples
ntrain = ceil(lengthOfExp*trainfrac);  
nRawTrain =  ceil(lengthOfExpRaw*trainfrac)

% number of test samples
ntest = lengthOfExp - ntrain; 
nRawtest = lengthOfExpRaw - nRawTrain;
 
% time indices for test
iitest = 1:ntest;
iiRawTest = 1:nRawtest;
 
% time indices for training
iitrain = ntest+1:lengthOfExp;  
iiRawTrain = nRawtest + 1:lengthOfExpRaw;
 
% training stimulus
stimtrain = scaledStimulus(iitrain);
% test stimulus
stimtest = scaledStimulus(iitest);

% Train spikes
spstrain = scaledSpikes(iitrain);
spsRawTrain = rawSpikesVector(iiRawTrain);

% Test spikes
spstest =  scaledSpikes(iitest);
spsRawTest = rawSpikesVector(iiRawTest);

spsCoupleddRawTrain = zeros(numOfCoupledNeurons, nRawTrain);
spsCoupleddRawTest = zeros(numOfCoupledNeurons, nRawtest);

for i = 1:numOfCoupledNeurons
    spsCoupleddRawTrain(i,:) = neuronRawStimulus(i,iiRawTrain);
    spsCoupleddRawTest(i,:) = neuronRawStimulus(i,iiRawTest);
end
% Print num of spikes 
fprintf('Taining: %d spikes', sum(spstrain));
fprintf('Testing: %d spikes', sum(spstest));
 
% Define the wanted size of stimulus filter before spike
filterSizeBeforeSpike = 200;
 
%% Post spike base vectors

% Define number of base vectors for post spike filter
numOfBaseVectors = 10;
 
% Define parameters for post spike base vectors
lastPeak = 1;
dt = 0.01;
hpeaks = [0.01 lastPeak];
b = 0.5;
 
% build post spike BaseVectors
[postSpiketimeVector,postSpikeBaseVectors, originalBaseVectors] = buildBaseVectorsForPostSpikeAndCoupling(numOfBaseVectors,dt,hpeaks, b);
 
% Update the size after base vectors build(Can be changed)
numOfBaseVectors = size(postSpikeBaseVectors,2);

% Plot base vectors
figure();
plot(postSpikeBaseVectors);
title('Base vectors for post spike history');
xlabel('Time after spike');

%% Design Matrix build

% We build the stimulus design matrix for train data
trainStimulusDesignMatrix = buildStimulusDesignMatrix(filterSizeBeforeSpike, stimtrain);

% We build the stimulus design matrix for test data
testStimulusDesignMatrix = buildStimulusDesignMatrix(filterSizeBeforeSpike, stimtest);

% We calculate the cell STA(With no zero-mean)
cellSTA = calculateSTA(trainStimulusDesignMatrix,spstrain);

% We build spike history design matrix
trainSpikeHistoryDesignMatrix = buildSpikeHistoryDesignMatrix(numOfBaseVectors,numOfCoupledNeurons, postSpikeBaseVectors, length(spstrain), spsRawTrain, wantedSampFactor, spsCoupleddRawTrain);
testSpikeHistoryDesignMatrix = buildSpikeHistoryDesignMatrix(numOfBaseVectors,numOfCoupledNeurons, postSpikeBaseVectors,  length(spstest), spsRawTest, wantedSampFactor, spsCoupleddRawTest);

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
cellPostSpike = zeros(1,numOfBaseVectors);
learnedParameters = [cellSTA' cellPostSpike];

% This matrix computes differences between adjacent coeffs
Dx1 = spdiags(ones(filterSizeBeforeSpike - 1,1)*[-1 1],0:1,filterSizeBeforeSpike-2,filterSizeBeforeSpike - 1);

% computes squared diffs(in order to get (w1 - w2).^2 penalty
Dx = Dx1'*Dx1; 

% Select lambda smoothing penalty by cross-validation 
 % grid of lambda values (ridge parameters)
lambdavals = 2.^(1:20);
nlambda = length(lambdavals);

% Embed Dx matrix in matrix with one extra row/column for constant coeff
D = blkdiag(0,Dx); 

% Allocate space for train and test errors
negLtrain = zeros(nlambda,1); 
negLogTest = zeros(nlambda,1);  
lambdaLearrnedParameters = zeros(length(learnedParameters),nlambda);


%% Run optimization problem with diffrent lambdas

figure();
hold on;

% Run for each lambda, and learn the the parameters
for i = 1:nlambda
    
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
    negLtrain(i) = Loss_GLM_logli_exp(learnedParameters, dataForLearnning);
    negLogTest(i) = Loss_GLM_logli_exp(learnedParameters, dataForTesting);
    
    % We plot the learned STA of current iteration
    plot(learnedParameters(1:filterSizeBeforeSpike));
    xlabel('time before spike');drawnow; pause(.5);
    ylabel('intensity');
end

hold off;

% Get the minimum log likelihood index
[~,imin] = min(negLogTest);

% Calculate the spike history vector based on the parameters learned
spikeHistoryVector = lambdaLearrnedParameters(end - numOfBaseVectors + 1 :end,imin)' * postSpikeBaseVectors';


%% Plot learned estimators

figure();

% Plot STA estimator
subplot(4,1,1);
plot(lambdaLearrnedParameters(1:filterSizeBeforeSpike,imin));
hold on;
plot(cellSTA(1:end) - mean(cellSTA));
legend('learned STA','Expiriment STA');
xlabel('Time before spike');
ylabel('intensity');
title('STA estimatror');

% Plot leaned spike history filter
subplot(4,1,2);
plot(spikeHistoryVector);
title('Spike history filter');
xlabel('Time after spike spike');
ylabel('factor to fire');

% Train likelihood
subplot(4,1,3);
plot(-negLtrain);
title('train likelihood');
xlabel('lambda factor');
ylabel('log likelihood');

% Test likelihood
subplot(4,1,4);
plot(-negLogTest);
title('test likelihood');
xlabel('lambda factor');
ylabel('log likelihood');
end
