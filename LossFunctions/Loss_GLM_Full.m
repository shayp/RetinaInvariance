function [logli, dL, H] = Loss_GLM_Full(learnedParameters,dataForLearnning)
% [neglogli, dL, H] = Loss_GLM_logli_exp(prs)
%(Taken from pillowLab)
% Compute negative log-likelihood of data undr the GLM model with
% exponential nonlinearity (with standard linear parametrization of stim filter) 

% Extract some vals from dataForLearnning:
stimulusFilterSize = dataForLearnning.stimulusFiltetSize;
% absolute bin size for spike train (in sec)
binSizeInSecond = dataForLearnning.binSizeInSecond;

% Unpack GLM prs;
stimulusFilter = learnedParameters(1:stimulusFilterSize);
meanFiringRate = learnedParameters(stimulusFilterSize + 1);
postspikehistoryFilters = learnedParameters(stimulusFilterSize + 2:end);
interpMatrix = dataForLearnning.interpMatrix;

spikeIndexes = dataForLearnning.spikeIndexes;
nsp = sum(dataForLearnning.spikesTrain);

% Extract some other stuff we'll use a lot
stimulusDesignMatrix = dataForLearnning.stimulusDesignMatrix; % stimulus design matrix
spikeHistoryDesignMatrix = dataForLearnning.spikeHistoryDesignMatrix';    % spike history design matrix
dataLen = dataForLearnning.dataLen;   % number of bins in spike train vector

% -------- Compute sum of filter reponses -----------------------
 linearFilter = interpMatrix * (stimulusDesignMatrix * stimulusFilter) + spikeHistoryDesignMatrix * postspikehistoryFilters + meanFiringRate; 

 % ---------  Compute output of nonlinearity  ------------------------
expValue = exp(linearFilter);

% ---------  Compute log-likelihood ---------------------------------
Trm0 = sum(expValue)* binSizeInSecond;  % non-spike term
Trm1 =  -sum(linearFilter(spikeIndexes)); % spike term
logli = Trm0 + Trm1;

% ---------  Compute Gradient -----------------
    % Non-spiking terms (Term 1)
    dLdStimulusFilter0 = (expValue' * interpMatrix * stimulusDesignMatrix)';
   
    dLdMeanFiringRate0 = sum(expValue);
    
    dLdSpikeHistoryFilter0 = (expValue' * spikeHistoryDesignMatrix)';
    
    Msp = interpMatrix(spikeIndexes,:); % interpolation matrix just for spike bins

    % Spiking terms (Term 2)
    dLdStimulusFilter1 = (sum(Msp*stimulusDesignMatrix))';
    
    dLdMeanFiringRate1 = nsp;
    dLdSpikeHistoryFilter1 = sum(spikeHistoryDesignMatrix(spikeIndexes,:),1)';

    % Combine terms
    dLdStimulusFilter = dLdStimulusFilter0 * binSizeInSecond  - dLdStimulusFilter1;
    dLdMeanFiringRate = dLdMeanFiringRate0*binSizeInSecond - dLdMeanFiringRate1;
    dLdSpikeHistoryFilter = dLdSpikeHistoryFilter0 * binSizeInSecond - dLdSpikeHistoryFilter1;
    
    dL = [dLdStimulusFilter; dLdMeanFiringRate; dLdSpikeHistoryFilter];
    
    rrdiag = spdiags(expValue, 0, dataLen, dataLen);
    hInterp = rrdiag * interpMatrix;
    
    Hk = (stimulusDesignMatrix' * (interpMatrix' * hInterp) * stimulusDesignMatrix) * binSizeInSecond;
    Hb = dLdMeanFiringRate0 * binSizeInSecond;
    Hh = (spikeHistoryDesignMatrix' * (bsxfun(@times,spikeHistoryDesignMatrix,expValue))) * binSizeInSecond;
    
    Hkb = (sum(hInterp,1) * stimulusDesignMatrix)' * binSizeInSecond;
    Hkh = ((spikeHistoryDesignMatrix'*hInterp)*stimulusDesignMatrix * binSizeInSecond)';
    Hhb = (expValue' * spikeHistoryDesignMatrix)' * binSizeInSecond;
    H = [[Hk Hkb Hkh]; [Hkb' Hb Hhb']; [Hkh' Hhb Hh]];

end
