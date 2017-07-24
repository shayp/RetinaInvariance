function [result_GLM_Full, result_GLM_Partial, result_LN] = ...
    runLearningModels(trainStimulusDesignMatrix, testStimulusDesignMatrix,...
    trainSpikeHistoryDesignMatrix, testSpikeHistoryDesignMatrix,...
    interpMatrixTrain, interpMatrixTest, trainSpikesTrain, testSpikesTrain, stimulusFilterParamsSize, binsInSecond,initStimulusFilter,...
    numOfCoupledNeurons, numOfBaseVectors, postSpikeBaseVectors, stimulusFilterSizeForSimulation)
%% Initlaization
addpath('./LossFunctions', './preProcessTools');
load('globalParams');
%% Optimization problem params init

% We create struct in order to gather data for the optimization
% problem(tesing and training)
dataForLearnning.stimulusFiltetSize = stimulusFilterParamsSize;
dataForLearnning.binSizeInSecond = 1 / binsInSecond ;
dataForLearnning.stimulusDesignMatrix = trainStimulusDesignMatrix;
dataForLearnning.spikeHistoryDesignMatrix = trainSpikeHistoryDesignMatrix;
dataForLearnning.dataLen = length(trainSpikesTrain);
dataForLearnning.spikesTrain = trainSpikesTrain;
dataForLearnning.interpMatrix = interpMatrixTrain;

dataForTesting.stimulusFiltetSize = stimulusFilterParamsSize;
dataForTesting.binSizeInSecond = 1 / binsInSecond ;
dataForTesting.stimulusDesignMatrix = testStimulusDesignMatrix;
dataForTesting.spikeHistoryDesignMatrix = testSpikeHistoryDesignMatrix;
dataForTesting.dataLen = length(testSpikesTrain);
dataForTesting.spikesTrain = testSpikesTrain;
dataForTesting.interpMatrix = interpMatrixTest;

% set options for optimizatiom problem
opts = optimset('Gradobj','on','Hessian','on','display','off');

% Set start parameters for the optimization problem
cellPostSpike =  0.0001 * ones(size(trainSpikeHistoryDesignMatrix,1), 1);
meanFiringRate = 0.0001;

% This matrix computes differences between adjacent coeffs
Dx1 = spdiags(ones(stimulusFilterParamsSize,1)*[-1 1],0:1,stimulusFilterParamsSize - 1,stimulusFilterParamsSize);
% computes squared diffs(in order to get (w1 - w2).^2 penalty
Dx = Dx1'*Dx1; 
% Select lambda smoothing penalty by cross-validation 
 % grid of lambda values (ridge parameters)
lambdavals = (2).^(6:17);
nlambda = length(lambdavals);


%% Run optimization problem with diffrent lambdas
params_LN_Bias = [initStimulusFilter' meanFiringRate]';
params_GLM_Full = [initStimulusFilter' cellPostSpike' meanFiringRate]';
params_GLM_Partial = [cellPostSpike' meanFiringRate]';
params_LN_No_Bias = initStimulusFilter;

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

timeSeries = linspace(-stimulusFilterSizeForSimulation * deltaT, 0, stimulusFilterParamsSize);
figure();
% Run for each lambda, and learn the the parameters
for i = 1:nlambda
    currentIteration = i
    % Compute the inverce covariance matrix with the current lambda
    Cinv = lambdavals(i) * Dx; % set inverse prior covariance
    
    % Learn  LN model with bias
    negLikelihhod = @(prs)Loss_LN_Bias(prs,dataForLearnning);
    lossfun = @(prs)Loss_posterior_LN_bias(prs,negLikelihhod,Cinv, stimulusFilterParamsSize);
    learned_LN_Bias(:,i) = fminunc(lossfun,params_LN_Bias,opts);
    LL_LN_Train(i) = Loss_LN_Bias(learned_LN_Bias(:,i), dataForLearnning);
    LL_LN_Test(i) = Loss_LN_Bias(learned_LN_Bias(:,i), dataForTesting);
    
    % Learn Full GLM
    negLikelihhod = @(prs)Loss_GLM_Full(prs,dataForLearnning);
    lossfun = @(prs)Loss_posterior_GLM(prs,negLikelihhod,Cinv, stimulusFilterParamsSize);
    learned_GLM_Full(:,i) = fminunc(lossfun,params_GLM_Full,opts);
    LL_GLM_Full_Train(i) = Loss_GLM_Full(learned_GLM_Full(:,i), dataForLearnning);
    LL_GLM_Full_Test(i) = Loss_GLM_Full(learned_GLM_Full(:,i), dataForTesting);
    
    % Learn LN no bias
    negLikelihhod = @(prs)Loss_LN__NoBias(prs,dataForLearnning);
    lossfun = @(prs)Loss_posterior_LN_NoBias(prs,negLikelihhod,Cinv);
    learned_LN_No_Bias(:,i) = fminunc(lossfun,params_LN_No_Bias,opts);
    
    % Learn GLM partial - using LN_No_Bias stimulus filter
    dataForLearnning.stimulusProjection = interpMatrixTrain * (trainStimulusDesignMatrix *  learned_LN_No_Bias(:,i));
    dataForTesting.stimulusProjection = interpMatrixTest * (testStimulusDesignMatrix *  learned_LN_No_Bias(:,i));
    lossfun = @(prs)Loss_GLM_Partial_History(prs,dataForLearnning);
    current_params_GLM_Partial = fminunc(lossfun,params_GLM_Partial,opts);
    learned_GLM_Partial(:,i) = [learned_LN_No_Bias(:,i)' current_params_GLM_Partial']';
    LL_GLM_Partial_Train(i) = Loss_GLM_Partial_History(current_params_GLM_Partial, dataForLearnning);
    LL_GLM_Partial_Test(i) = Loss_GLM_Partial_History(current_params_GLM_Partial, dataForTesting);

    clf;
    hold on;
    plot(timeSeries, initStimulusFilter,...
         timeSeries, learned_LN_Bias(1:stimulusFilterParamsSize,i),...
         timeSeries, learned_GLM_Full(1:stimulusFilterParamsSize,i),...
         timeSeries, learned_LN_No_Bias(1:stimulusFilterParamsSize,i));
    legend('STA', 'LN', 'GLM Full', 'GLM Partial');
    xlabel('Time before spike(s)');
    ylabel('intensity');
    title(['Stimulus filter estimator - iteration = :'  num2str(i)]);
    drawnow;
    hold off;
end

% Get the minimum log likelihood index
[~,imin_LN_Bias] = min(LL_LN_Test);
[~,imin_GLM_Full] = min(LL_GLM_Full_Test);
[~,imin_GLM_Partial] = min(LL_GLM_Partial_Test);
imin_LN_Bias
imin_GLM_Full
imin_GLM_Partial
bestParams_LN = learned_LN_Bias(:, imin_LN_Bias);
bestParams_GLM_Full = learned_GLM_Full(:, imin_GLM_Full);
bestParams_GLM_Partial = learned_GLM_Partial(:, imin_GLM_Partial);

for Index = 1:numOfCoupledNeurons
    result_GLM_Full.couplingFilters(Index,:) = bestParams_GLM_Full(stimulusFilterParamsSize + (Index  - 1) * numOfBaseVectors + 1 :stimulusFilterParamsSize + (Index) * numOfBaseVectors)' * postSpikeBaseVectors';
    result_GLM_Partial.couplingFilters(Index,:) = bestParams_GLM_Partial(stimulusFilterParamsSize + (Index  - 1) * numOfBaseVectors + 1 :stimulusFilterParamsSize + (Index) * numOfBaseVectors)' * postSpikeBaseVectors';
end

fineTimeScale = linspace(-deltaT * stimulusFilterSizeForSimulation, 0, stimulusFilterSizeForSimulation);
coarseTimeScale = linspace(-deltaT * stimulusFilterSizeForSimulation, 0, stimulusFilterParamsSize);

result_GLM_Full.StimulusFilter = interp1(coarseTimeScale, bestParams_GLM_Full(1:stimulusFilterParamsSize), fineTimeScale, 'spline');
result_GLM_Full.meanFiringRate = bestParams_GLM_Full(end);

result_GLM_Partial.StimulusFilter = interp1(coarseTimeScale, bestParams_GLM_Partial(1:stimulusFilterParamsSize), fineTimeScale, 'spline');
result_GLM_Partial.meanFiringRate = bestParams_GLM_Partial(end);

result_LN.StimulusFilter = interp1(coarseTimeScale, bestParams_LN(1:stimulusFilterParamsSize), fineTimeScale, 'spline');
result_LN.meanFiringRate = bestParams_LN(end);
% figure();
% plot(timeSeries, initStimulusFilter,...
%      timeSeries, bestParams_LN(1:stimulusFilterParamsSize),...
%      timeSeries, bestParams_GLM_Full(1:stimulusFilterParamsSize),...
%      timeSeries, bestParams_GLM_Partial(1:stimulusFilterParamsSize));
% legend('STA', 'LN', 'GLM Full', 'GLM Partial');
% xlabel('Time before spike(s)');
% ylabel('intensity');
% title('Choosen stimulus filter(before resize):');
% drawnow;

end
