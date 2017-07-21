function [logli, dL, H] = Loss_LN__NoBias(learnedParameters,dataForLearnning)
%(Taken from pillowLab)
% Compute negative log-likelihood of data undr the LN model with
% exponential nonlinearity (with standard linear parametrization of stim filter) 

% absolute bin size for spike train (in sec)
binSizeInSecond = dataForLearnning.binSizeInSecond;    

interpMatrix = dataForLearnning.interpMatrix;

% Unpack LN params;
stimulusFilter = learnedParameters;
spikesTrain = dataForLearnning.spikesTrain;
dataLen = dataForLearnning.dataLen;   % number of bins in spike train vector

% Extract some other stuff we'll use a lot
stimulusDesignMatrix = dataForLearnning.stimulusDesignMatrix; % stimulus design matrix

% -------- Compute sum of filter reponses -----------------------
 linearFilter = interpMatrix * (stimulusDesignMatrix * stimulusFilter);
% ---------  Compute output of nonlinearity  ------------------------
expValue = exp(linearFilter) * binSizeInSecond;

% ---------  Compute log-likelihood ---------------------------------
Trm0 = sum(expValue);  % non-spike term
Trm1 = -linearFilter' * spikesTrain; % spike term
logli = Trm0 + Trm1;

% ---------  Compute Gradient -----------------
if (nargout > 1)
    
    % Non-spiking terms (Term 1)
    dLdStimulusFilter0 = (expValue' * interpMatrix * stimulusDesignMatrix)';
    
    % Spiking terms (Term 2)
    dLdStimulusFilter1 = (interpMatrix * stimulusDesignMatrix)' * spikesTrain;
    
    % Combine terms
    dLdStimulusFilter = dLdStimulusFilter0  - dLdStimulusFilter1;
    
    dL = dLdStimulusFilter';
end

 if (nargout > 2)
    rrdiag = spdiags(expValue, 0, dataLen, dataLen); 
    hInterp = rrdiag * interpMatrix;
    H = (stimulusDesignMatrix' * (interpMatrix' * hInterp) * stimulusDesignMatrix) * binSizeInSecond;
 end
 
end
