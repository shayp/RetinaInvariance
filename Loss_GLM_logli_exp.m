function [logli, dL, H, Hk, Hkb, Hkh, Hb, Hhb, Hh] = Loss_GLM_logli_exp(learnedParameters,dataForLearnning)
% [neglogli, dL, H] = Loss_GLM_logli_exp(prs)
%(Taken from pillowLab)
% Compute negative log-likelihood of data undr the GLM model with
% exponential nonlinearity (with standard linear parametrization of stim filter) 
%
% Inputs:
%   prs = [kprs - weights for stimulus kernel
%          dc   - dc current injection
%          ihprs - weights on post-spike current];
% Outputs:
%      logli = negative log likelihood of spike train
%      dL = gradient with respect to prs
%      H = hessian

% Extract some vals from Xstruct (Opt Prs);
stimulusFilterSize = dataForLearnning.stimulusFiltetSize;   % total # params for k
binSizeInSecond = dataForLearnning.binSizeInSecond;           % absolute bin size for spike train (in sec)

% Unpack GLM prs;
stimulusFilter = learnedParameters(1:stimulusFilterSize);
postspikehistoryFilters = learnedParameters(stimulusFilterSize+1:end - 1);
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
    
    % Non-spiking terms (Term 1)
    dLdStimulusFilter0 = stimulusDesignMatrix' * expValue;
   
    dLdMeanFiringRate0 = sum(expValue);
    
    dLdSpikeHistoryFilter0 = spikeHistoryDesignMatrix' * expValue;
    
    % Spiking terms (Term 2)
    dLdStimulusFilter1 = stimulusDesignMatrix' * spikesTrain;
    
    dLdMeanFiringRate1 = nsp;
    dLdSpikeHistoryFilter1 = spikeHistoryDesignMatrix' * spikesTrain;
    
    % Combine terms
    dLdStimulusFilter = dLdStimulusFilter0 * binSizeInSecond  - dLdStimulusFilter1;
    dLdMeanFiringRate = dLdMeanFiringRate0*binSizeInSecond - dLdMeanFiringRate1;
    dLdSpikeHistoryFilter = dLdSpikeHistoryFilter0 * binSizeInSecond - dLdSpikeHistoryFilter1;
    
    dL = [dLdStimulusFilter' dLdSpikeHistoryFilter' dLdMeanFiringRate];
end
 if (nargout > 2)
    ddrrdiag = spdiags(expValue, 0, dataLen, dataLen); 
    Hk = stimulusDesignMatrix' * bsxfun(@times,stimulusDesignMatrix,expValue) * binSizeInSecond;
    Hb = dLdMeanFiringRate0 * binSizeInSecond;
    Hh = spikeHistoryDesignMatrix' * bsxfun(@times,spikeHistoryDesignMatrix,expValue) * binSizeInSecond;  % Hh (h filter)
    Hkh = ((spikeHistoryDesignMatrix' * ddrrdiag) * stimulusDesignMatrix)' * binSizeInSecond;         % Hhk (cross-term)
    Hhb = (expValue' * spikeHistoryDesignMatrix)' * binSizeInSecond;
    Hkb = (sum(ddrrdiag,1) * stimulusDesignMatrix)' * binSizeInSecond;
    H = [[Hk Hkb Hkh]; [Hkb' Hb Hhb']; [Hkh' Hhb Hh]];

 end
end
