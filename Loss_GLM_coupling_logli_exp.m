function [logli, dL, H] = Loss_GLM_coupling_logli_exp(learnedParameters,dataForLearnning)
%(Taken from pillowLab)
% Compute negative log-likelihood of data undr the GLM model with
% exponential nonlinearity (with standard linear parametrization of stim filter) 
%

stimulusFilterSize = dataForLearnning.stimulusFiltetSize;   % total # params for k
binSizeInSecond = dataForLearnning.binSizeInSecond;           % absolute bin size for spike train (in sec)

% Unpack GLM prs;
stimulusFilter = dataForLearnning.stimulusFilter;
postspikehistoryFilters = learnedParameters(1:end - 1);
meanFiringRate = learnedParameters(end);
spikesTrain = dataForLearnning.spikesTrain;
nsp = sum(dataForLearnning.spikesTrain);
% Extract some other stuff we'll use a lot
stimulusDesignMatrix = dataForLearnning.stimulusDesignMatrix; % stimulus design matrix
spikeHistoryDesignMatrix = dataForLearnning.spikeHistoryDesignMatrix;    % spike history design matrix
dataLen = dataForLearnning.dataLen;   % number of bins in spike train vector
% -------- Compute sum of filter reponses -----------------------
 linearFilter = stimulusDesignMatrix*stimulusFilter' + spikeHistoryDesignMatrix * postspikehistoryFilters'; 
 linearFilter = linearFilter + meanFiringRate;

 % ---------  Compute output of nonlinearity  ------------------------
expValue = exp(linearFilter);

% ---------  Compute log-likelihood ---------------------------------
Trm0 = sum(expValue)* binSizeInSecond;  % non-spike term
Trm1 =  -linearFilter' * spikesTrain; % spike term
logli = Trm0 + Trm1;

% ---------  Compute Gradient -----------------
if (nargout > 1)
    
    dLdMeanFiringRate0 = sum(expValue);
    
    dLdSpikeHistoryFilter0 = spikeHistoryDesignMatrix' * expValue;
    
    dLdMeanFiringRate1 = nsp;
    dLdSpikeHistoryFilter1 = spikeHistoryDesignMatrix' * spikesTrain;
    
    % Combine terms
    dLdMeanFiringRate = dLdMeanFiringRate0*binSizeInSecond - dLdMeanFiringRate1;
    dLdSpikeHistoryFilter = dLdSpikeHistoryFilter0 * binSizeInSecond - dLdSpikeHistoryFilter1;
    
    dL = [dLdSpikeHistoryFilter' dLdMeanFiringRate];
end
 if (nargout > 2)
    ddrrdiag = spdiags(expValue, 0, dataLen, dataLen); 
    Hb = dLdMeanFiringRate0 * binSizeInSecond;
    Hh = spikeHistoryDesignMatrix' * bsxfun(@times,spikeHistoryDesignMatrix,expValue) * binSizeInSecond;  % Hh (h filter)
    Hhb = (expValue' * spikeHistoryDesignMatrix)' * binSizeInSecond;
    H = [[Hh Hhb]; [Hhb' Hb]];

 end
end
