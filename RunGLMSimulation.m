function response = RunGLMSimulation(numOfNeurons, Stimulus, Filters, stimulusFilterLength, couplingFilterLength,deltaT)

maxFilterLength = max(stimulusFilterLength, couplingFilterLength);
simulationLength = length(Stimulus);
response = zeros(numOfNeurons, simulationLength + maxFilterLength);
projectedStimulus = zeros(numOfNeurons, simulationLength + maxFilterLength);
nCount = 0;
nbinsPerEval = 100;  % Default number of bins to update for each spike

for neuronIndex = 1:numOfNeurons
%     if Filters(neuronIndex).meanFiringRate < 0
%         Filters(neuronIndex).meanFiringRate = 0;
%     end
    projectedStimulus(neuronIndex, stimulusFilterLength + 1:end) = conv(Stimulus, Filters(neuronIndex).StimulusFilter, 'same');
    projectedStimulus(neuronIndex,:) = projectedStimulus(neuronIndex,:) + Filters(neuronIndex).meanFiringRate;
end

rprev = 0;
nsp = 0;
tspnext = exprnd(1);  % time of next spike (in rescaled time)

currentBin= maxFilterLength + 1;
% while currentBin <= simulationLength
%     iinxt = currentBin:min(currentBin+nbinsPerEval-1,simulationLength);
%     rrnxt = exp(projectedStimulus(1,iinxt)) * deltaT; % Cond Intensity
%     rrcum = cumsum(rrnxt)+ rprev; % integrated cond intensity
%     if (tspnext >= rrcum(end)) % No spike in this window
%         currentBin = iinxt(end)+1;
%         rprev = rrcum(end);
%     else   % Spike!
%         ispk = iinxt(find(rrcum>=tspnext, 1, 'first')); % time bin where spike occurred
%         nsp = nsp+1;
%         response(1,ispk) = 1; % spike time
%          mxi = min(simulationLength, ispk+couplingFilterLength); % max time affected by post-spike kernel
%          iiPostSpk = ispk+1:mxi; % time bins affected by post-spike kernel
%         projectedStimulus(1,iiPostSpk) = projectedStimulus(1,iiPostSpk)+  Filters(1).couplingFilters(1,mxi-ispk);
% %         end
%         tspnext = exprnd(1);  % draw next spike time
%         rprev = 0; % reset integrated intensity
%         currentBin = ispk+1;  % Move to next bin
%         % --  Update # of samples per iter ---
%         muISI = currentBin/nsp;
%         nbinsPerEval = max(20, round(1.5*muISI)); 
%     end
% end
for i = maxFilterLength + 1:simulationLength
    for neuronIndex = 1:numOfNeurons
        projectionTrm = projectedStimulus(neuronIndex,i);
        for couplingIndex = 1:numOfNeurons
            projectionTrm = projectionTrm + Filters(neuronIndex).couplingFilters(couplingIndex,:) * response(couplingIndex, i - couplingFilterLength:i - 1)';
        end
        curentLambda = exp(projectionTrm) * deltaT;
        sample = poissrnd(curentLambda);
        if sample > 1
            nCount  = nCount + 1;
            sample = 1;
        end
      response(neuronIndex, i) = sample;
        
    end
end
 response = response(:,maxFilterLength + 1:end);

end