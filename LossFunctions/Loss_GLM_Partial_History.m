function [logli, dL, H] = Loss_GLM_Partial_History(learnedParameters,dataForLearnning)
%(Taken from pillowLab)
% Compute negative log-likelihood of data undr the GLM model with
% exponential nonlinearity (with standard linear parametrization of stim filter) 
%

% absolute bin size for spike train (in sec)
binSizeInSecond = dataForLearnning.binSizeInSecond;           

% Unpack GLM prs;
constantStimulusProjection = dataForLearnning.stimulusProjection;
postspikehistoryFilters = learnedParameters(1:end);
spikesTrain = dataForLearnning.spikesTrain;

% Extract some other stuff we'll use a lot
spikeHistoryDesignMatrix = dataForLearnning.spikeHistoryDesignMatrix';    % spike history design matrix

% -------- Compute sum of filter reponses -----------------------
linearFilter = constantStimulusProjection + spikeHistoryDesignMatrix * postspikehistoryFilters; 

 % ---------  Compute output of nonlinearity  ------------------------
expValue = exp(linearFilter);

% ---------  Compute log-likelihood ---------------------------------
Trm0 = sum(expValue) * binSizeInSecond;  % non-spike term
Trm1 =  -linearFilter' * spikesTrain; % spike term
logli = Trm0 + Trm1;

% ---------  Compute Gradient -----------------
if (nargout > 1)
    dLdSpikeHistoryFilter0 = (spikeHistoryDesignMatrix' * expValue);
    dLdSpikeHistoryFilter1 = spikeHistoryDesignMatrix' * spikesTrain;
    
    % Combine terms
    dLdSpikeHistoryFilter = dLdSpikeHistoryFilter0 * binSizeInSecond - dLdSpikeHistoryFilter1;
    
    dL = dLdSpikeHistoryFilter';
end
 if (nargout > 2)
    H = spikeHistoryDesignMatrix' * bsxfun(@times,spikeHistoryDesignMatrix,expValue) * binSizeInSecond;  % Hh (h filter)
 end
end
