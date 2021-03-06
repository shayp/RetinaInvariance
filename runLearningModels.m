function [result_GLM_Full, result_LN] = ...
    runLearningModels(trainStimulusDesignMatrix, testStimulusDesignMatrix,...
    trainSpikeHistoryDesignMatrix, testSpikeHistoryDesignMatrix,...
    interpMatrixTrain, interpMatrixTest, trainSpikesTrain, testSpikesTrain, stimulusFilterParamsSize, binsInSecond,initStimulusFilter,...
    numOfCoupledNeurons, numOfSpikeHistoryBaseVectors, numOfCouplingBaseVectors, spikeHistoryBaseVectors, couplingBaseVectors, stimulusFilterSizeForSimulation)
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
dataForLearnning. spikeIndexes = find(trainSpikesTrain);

dataForTesting.stimulusFiltetSize = stimulusFilterParamsSize;
dataForTesting.binSizeInSecond = 1 / binsInSecond ;
dataForTesting.stimulusDesignMatrix = testStimulusDesignMatrix;
dataForTesting.spikeHistoryDesignMatrix = testSpikeHistoryDesignMatrix;
dataForTesting.dataLen = length(testSpikesTrain);
dataForTesting.spikesTrain = testSpikesTrain;
dataForTesting.interpMatrix = interpMatrixTest;
dataForTesting. spikeIndexes = find(testSpikesTrain);

% set options for optimizatiom problem
%opts = optimset('Gradobj','on','Hessian','on','display','off');
opts = optimset('Gradobj','on','Hessian','on');

% Set start parameters for the optimization problem
cellPostSpike =  zeros(size(trainSpikeHistoryDesignMatrix,2), 1);
cellPostSpike(1:end) = linspace(-10, 0, length(cellPostSpike))';
meanFiringRate = 0;

% This matrix computes differences between adjacent coeffs
Dx1 = spdiags(ones(stimulusFilterParamsSize,1)*[-1 1],0:1,stimulusFilterParamsSize - 1,stimulusFilterParamsSize);

% computes squared diffs(in order to get (w1 - w2).^2 penalty
Dx = Dx1'*Dx1; 

% Select lambda smoothing penalty by cross-validation 
 % grid of lambda values (ridge parameters)
lambdavals = (2).^(7:15);

nlambda = length(lambdavals);

%% Run optimization problem with diffrent lambdas
params_LN_Bias = [initStimulusFilter; meanFiringRate];
params_GLM_Full = [initStimulusFilter; meanFiringRate; cellPostSpike];

learned_LN_Bias = zeros(length(params_LN_Bias),nlambda);
learned_GLM_Full = zeros(length(params_GLM_Full),nlambda);

% Allocate space for train and test errors
LL_GLM_Full_Train = zeros(nlambda,1);
LL_GLM_Full_Test = zeros(nlambda,1); 
LL_LN_Train = zeros(nlambda,1);
LL_LN_Test = zeros(nlambda,1);

timeSeries = linspace(-stimulusFilterSizeForSimulation * deltaT, 0, stimulusFilterParamsSize);
figure();
% Run for each lambda, and learn the the parameters
for i = 1:nlambda
    currentLambda = lambdavals(i)
    % Compute the inverce covariance matrix with the current lambda
    Cinv = lambdavals(i) * Dx; % set inverse prior covariance
    
    disp('************ LN Bias Learning *****************');
    % Learn  LN model with bias
    negLikelihhod = @(prs)Loss_LN_Bias(prs,dataForLearnning);
    lossfun = @(prs)Loss_posterior_LN_bias(prs,negLikelihhod,Cinv, stimulusFilterParamsSize);
    learned_LN_Bias(:,i) = fminunc(lossfun,params_LN_Bias,opts);
    LL_LN_Train(i) = Loss_LN_Bias(learned_LN_Bias(:,i), dataForLearnning);
    LL_LN_Test(i) = Loss_LN_Bias(learned_LN_Bias(:,i), dataForTesting);

    
    disp('************ GLM Full Learning *****************');
    % Learn Full GLM
    negLikelihhod = @(GlmParams)Loss_GLM_Full(GlmParams,dataForLearnning);
    lossfun = @(postGlmParams)Loss_posterior_GLM(postGlmParams,negLikelihhod,Cinv, stimulusFilterParamsSize);
    learned_GLM_Full(:,i) = fminunc(lossfun,params_GLM_Full,opts);
    LL_GLM_Full_Train(i) = Loss_GLM_Full(learned_GLM_Full(:,i), dataForLearnning);
    LL_GLM_Full_Test(i) = Loss_GLM_Full(learned_GLM_Full(:,i), dataForTesting);
    
    clf;
    hold on;
    subplot(1,2,1);
    plot(exp(learned_GLM_Full(stimulusFilterParamsSize + 2:stimulusFilterParamsSize + numOfSpikeHistoryBaseVectors + 1, i)' * spikeHistoryBaseVectors'));
    legend('Coupling Filter');
    xlabel('Time after spike spike');
    ylabel('intensity');
    title(['History filter estimator']);
    drawnow;
    hold off;
    hold on;
    subplot(1,2,2);
    plot(timeSeries,initStimulusFilter,...
         timeSeries, learned_GLM_Full(1:stimulusFilterParamsSize, i));
    legend('STA', 'GLM');
    xlabel('Time before  spike(s)');
    ylabel('intensity');
    title(['Stimulus filter estimator - iteration = :'  num2str(i)]);
    drawnow;
    hold off;
end

% Get the minimum log likelihood index
[~,imin_LN_Bias] = min(LL_LN_Test);
[~,imin_GLM_Full] = min(LL_GLM_Full_Test);
imin_LN_Bias
imin_GLM_Full
% imin_GLM_Partial
bestParams_LN = learned_LN_Bias(:, imin_LN_Bias);
bestParams_GLM_Full = learned_GLM_Full(:, imin_GLM_Full);
result_GLM_Full.spikeHistoryFilter = bestParams_GLM_Full(stimulusFilterParamsSize + 2 :stimulusFilterParamsSize + numOfSpikeHistoryBaseVectors + 1)' * spikeHistoryBaseVectors';
for Index = 1:numOfCoupledNeurons
    result_GLM_Full.couplingFilters(Index,:) = bestParams_GLM_Full(stimulusFilterParamsSize + numOfSpikeHistoryBaseVectors + (Index  - 1) * numOfCouplingBaseVectors + 2 :stimulusFilterParamsSize + numOfSpikeHistoryBaseVectors + (Index) * numOfCouplingBaseVectors + 1)' * couplingBaseVectors';
end

fineTimeScale = linspace(-deltaT * stimulusFilterSizeForSimulation, 0, stimulusFilterSizeForSimulation);
coarseTimeScale = linspace(-deltaT * stimulusFilterSizeForSimulation, 0, stimulusFilterParamsSize);

result_GLM_Full.StimulusFilter = interp1(coarseTimeScale, bestParams_GLM_Full(1:stimulusFilterParamsSize), fineTimeScale, 'spline');
result_GLM_Full.meanFiringRate = bestParams_GLM_Full(stimulusFilterParamsSize + 1);

result_LN.StimulusFilter = interp1(coarseTimeScale, bestParams_LN(1:stimulusFilterParamsSize), fineTimeScale, 'spline');
result_LN.meanFiringRate = bestParams_LN(end);
end
