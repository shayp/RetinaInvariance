function [logli, dL] = Loss_GLM_logli_exp(learnedParameters,dataForLearnning)
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
firingRateParam = learnedParameters(stimulusFilterSize+1);
postspikehistoryFilters = learnedParameters(stimulusFilterSize+2:end);

% Extract some other stuff we'll use a lot
stimulusDesignMatrix = dataForLearnning.stimulusDesignMatrix; % stimulus design matrix
spikeHistoryDesignMatrix = dataForLearnning.spikeHistoryDesignMatrix;    % spike history design matrix
spikesoccurence = dataForLearnning.spikesccurence;   % binary spike vector
%M = Xstruct.Minterp;   % matrix for interpolating from stimulus bins to spike train bins
dataLen = dataForLearnning.dataLen;   % number of bins in spike train vector
nsp = sum(spikesoccurence);     % number of spikes

% -------- Compute sum of filter reponses -----------------------
linearFilter = stimulusDesignMatrix*stimulusFilter + spikeHistoryDesignMatrix*postspikehistoryFilters + firingRateParam; % stim-dependent + spikehist-dependent inputs

% ---------  Compute output of nonlinearity  ------------------------
expValue = exp(linearFilter);

% ---------  Compute log-likelihood ---------------------------------
Trm1 = sum(expValue)*binSizeInSecond;  % non-spike term
Trm2 = -sum(linearFilter(spikesoccurence)); % spike term
logli = Trm1 + Trm2;

% ---------  Compute Gradient -----------------
if (nargout > 1)
    
    % Non-spiking terms (Term 1)
    dLdStimulusFilter0 = (expValue'*stimulusDesignMatrix)';
    dLdMeanFiringRate0 = sum(expValue);
    dLdSpikeHistoryFilter0 = (expValue'*spikeHistoryDesignMatrix)'; 
    
    % Spiking terms (Term 2)
    dLdStimulusFilter1 = (sum(stimulusDesignMatrix))';
    dLdMeanFiringRate1 = nsp;
    dLdSpikeHistoryFilter1 = sum(spikeHistoryDesignMatrix(spikesoccurence,:),1)';

    % Combine terms
    dLdStimulusFilter = dLdStimulusFilter0 * binSizeInSecond - dLdStimulusFilter1;
    dLdMeanFiringRate = dLdMeanFiringRate0*binSizeInSecond - dLdMeanFiringRate1;
    dLdSpikeHistoryFilter = dLdSpikeHistoryFilter0*binSizeInSecond - dLdSpikeHistoryFilter1;
    
    dL = [dLdStimulusFilter; dLdMeanFiringRate; dLdSpikeHistoryFilter];  
end
