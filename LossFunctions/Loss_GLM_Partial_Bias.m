function [logli, dL, H] = Loss_GLM_Partial_Bias(learnedParameters,dataForLearnning)
%(Taken from pillowLab)
% Compute negative log-likelihood of data undr the GLM model with
% exponential nonlinearity (with standard linear parametrization of stim filter) 
%

% absolute bin size for spike train (in sec)
binSizeInSecond = dataForLearnning.binSizeInSecond;           

% Unpack GLM prs;
constantStimulusProjection = dataForLearnning.stimulusProjection;
constantPostSpike = dataForLearnning.postSpikeHistoryVector;
meanFiringRate = learnedParameters(1);
spikesTrain = dataForLearnning.spikesTrain;
nsp = sum(dataForLearnning.spikesTrain);

% -------- Compute sum of filter reponses -----------------------
linearFilter = constantStimulusProjection +constantPostSpike + meanFiringRate; 

 % ---------  Compute output of nonlinearity  ------------------------
expValue = exp(linearFilter);

% ---------  Compute log-likelihood ---------------------------------
Trm0 = sum(expValue) * binSizeInSecond;  % non-spike term
Trm1 =  -linearFilter' * spikesTrain; % spike term
logli = Trm0 + Trm1;

% ---------  Compute Gradient -----------------
if (nargout > 1)

    dLdMeanFiringRate0 = sum(expValue);
    dLdMeanFiringRate1 = nsp;

    dLdMeanFiringRate = dLdMeanFiringRate0*binSizeInSecond - dLdMeanFiringRate1;

    dL = dLdMeanFiringRate;
end
 if (nargout > 2)
    H = dLdMeanFiringRate0 * binSizeInSecond;
 end
end
