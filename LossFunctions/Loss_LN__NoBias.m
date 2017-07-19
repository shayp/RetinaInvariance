function [logli, dL, H] = Loss_LN__NoBias(learnedParameters,dataForLearnning)
%(Taken from pillowLab)
% Compute negative log-likelihood of data undr the LN model with
% exponential nonlinearity (with standard linear parametrization of stim filter) 

% absolute bin size for spike train (in sec)
binSizeInSecond = dataForLearnning.binSizeInSecond;    

% Unpack LN params;
stimulusFilter = learnedParameters;
spikesTrain = dataForLearnning.spikesTrain;

% Extract some other stuff we'll use a lot
stimulusDesignMatrix = dataForLearnning.stimulusDesignMatrix; % stimulus design matrix

% -------- Compute sum of filter reponses -----------------------
 linearFilter = stimulusDesignMatrix*stimulusFilter;
% ---------  Compute output of nonlinearity  ------------------------
expValue = exp(linearFilter) * binSizeInSecond;

% ---------  Compute log-likelihood ---------------------------------
Trm0 = sum(expValue);  % non-spike term
Trm1 = -linearFilter' * spikesTrain; % spike term
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
