function [logli, dL, H] = Loss_LN_logli_exp(learnedParameters,dataForLearnning)
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
spikesTrain = dataForLearnning.spikesTrain;

% Extract some other stuff we'll use a lot
stimulusDesignMatrix = dataForLearnning.stimulusDesignMatrix; % stimulus design matrix
spikesoccurence = dataForLearnning.spikesccurence;   % binary spike vector
dataLen = dataForLearnning.dataLen;   % number of bins in spike train vector
nsp = sum(spikesoccurence);     % number of spikes

% -------- Compute sum of filter reponses -----------------------
 linearFilter = stimulusDesignMatrix*stimulusFilter';
% ---------  Compute output of nonlinearity  ------------------------
expValue = exp(linearFilter) * binSizeInSecond;

% ---------  Compute log-likelihood ---------------------------------
Trm0 = sum(expValue);  % non-spike term
Trm1 = -1 * linearFilter' * spikesTrain; % spike term
logli = Trm0 + Trm1;

% ---------  Compute Gradient -----------------
if (nargout > 1)
    
    % Non-spiking terms (Term 1)
    dLdStimulusFilter0 = stimulusDesignMatrix' * expValue;
    
    % Spiking terms (Term 2)
    dLdStimulusFilter1 = stimulusDesignMatrix' * spikesTrain;
    
    % Combine terms
    dLdStimulusFilter = dLdStimulusFilter0  - dLdStimulusFilter1;
    
    dL = dLdStimulusFilter';
end
 if (nargout > 2)

   H = stimulusDesignMatrix' * bsxfun(@times,stimulusDesignMatrix,expValue);
 end
end