function [logli, dL, H, Hk, Hkb, Hb] = Loss_LN_Bias(learnedParameters,dataForLearnning)
%(Taken from pillowLab)
% Compute negative log-likelihood of data undr the LM model with
% exponential nonlinearity (with standard linear parametrization of stim filter) 

% Extract some vals from dataForLearnning:
stimulusFilterSize = dataForLearnning.stimulusFiltetSize;
% absolute bin size for spike train (in sec)
binSizeInSecond = dataForLearnning.binSizeInSecond;
interpMatrix = dataForLearnning.interpMatrix;
% Unpack GLM prs;
stimulusFilter = learnedParameters(1:stimulusFilterSize);
meanFiringRate = learnedParameters(end);

spikesTrain = dataForLearnning.spikesTrain;
nsp = sum(dataForLearnning.spikesTrain);
% Extract some other stuff we'll use a lot
stimulusDesignMatrix = dataForLearnning.stimulusDesignMatrix; % stimulus design matrix
dataLen = dataForLearnning.dataLen;   % number of bins in spike train vector

% -------- Compute sum of filter reponses -----------------------
 linearFilter = interpMatrix * (stimulusDesignMatrix * stimulusFilter); 
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
    dLdStimulusFilter0 = (expValue' * interpMatrix * stimulusDesignMatrix)';
   
    dLdMeanFiringRate0 = sum(expValue);
        
    % Spiking terms (Term 2)
    dLdStimulusFilter1 = (spikesTrain' * interpMatrix * stimulusDesignMatrix)';
    
    dLdMeanFiringRate1 = nsp;
    
    % Combine terms
    dLdStimulusFilter = dLdStimulusFilter0 * binSizeInSecond  - dLdStimulusFilter1;
    dLdMeanFiringRate = dLdMeanFiringRate0*binSizeInSecond - dLdMeanFiringRate1;
    
    dL = [dLdStimulusFilter' dLdMeanFiringRate];
end
 if (nargout > 2)
    rrdiag = spdiags(expValue, 0, dataLen, dataLen); 
    hInterp = rrdiag * interpMatrix;
    Hk = (stimulusDesignMatrix' * (interpMatrix' * hInterp) * stimulusDesignMatrix) * binSizeInSecond;
    Hb = dLdMeanFiringRate0 * binSizeInSecond;
    Hkb = (sum(hInterp,1) * stimulusDesignMatrix)' * binSizeInSecond;
    H = [[Hk Hkb]; [Hkb' Hb]];
 end
end