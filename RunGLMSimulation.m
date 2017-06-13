function response = RunGLMSimulation(numOfNeurons, Stimulus, Filters, stimulusFilterLength, couplingFilterLength,deltaT)
maxFilterLength = max(stimulusFilterLength, couplingFilterLength);
simulationLength = length(Stimulus);
response = zeros(numOfNeurons, simulationLength + maxFilterLength);
linearValue = zeros(numOfNeurons, simulationLength + maxFilterLength);
nCount = 0;
nbinsPerEval = 100;  % Default number of bins to update for each spike

for neuronIndex = 1:numOfNeurons
    linearValue(neuronIndex, stimulusFilterLength + 1:end) = conv(Stimulus, Filters(neuronIndex).StimulusFilter, 'same');
    linearValue(neuronIndex,:) = linearValue(neuronIndex,:) + Filters(neuronIndex).meanFiringRate ;
end

rprev = zeros(1,numOfNeurons);
nsp = zeros(1,numOfNeurons);
tspnext = exprnd(1,1,numOfNeurons);
currentBin= maxFilterLength + 1;
while currentBin <= simulationLength
    iinxt = currentBin:min(currentBin+nbinsPerEval-1,simulationLength);
    nii = length(iinxt);  % Number of bins
    rrnxt =  exp(linearValue(:,iinxt)) * deltaT; % Cond Intensity
    rrcum = cumsum(rrnxt'+[rprev;zeros(nii-1,numOfNeurons)],1);  % Cumulative intensity
    if all(tspnext >= rrcum(end,:)) % No spike in this window
            currentBin = iinxt(end)+1;
            rprev = rrcum(end,:);
    else % Spike!
        [ispks,jspks] =  find(rrcum>=repmat(tspnext,nii,1));
        spcells = unique(jspks(ispks == min(ispks))); % cell number(s)
        ispk = iinxt(min(ispks)); % time bin of spike(s)
        rprev = rrcum(min(ispks),:); % grab accumulated history to here
        % Record this spike
        mxi = min(maxFilterLength, ispk+couplingFilterLength); % determine bins for adding h current
        iiPostSpk = ispk+1:mxi;
         for ic = 1:length(spcells)
            icell = spcells(ic);
            nsp(icell) = nsp(icell)+1;
            response(icell, ispk) = 1;
            for i = 1:numOfNeurons
                linearValue(i,iiPostSpk) = linearValue(i,iiPostSpk) + Filters(i).couplingFilters(icell,1:mxi-ispk);
            end
            rprev(icell) = 0;  % reset this cell's integral
            tspnext(icell) = exprnd(1); % draw RV for next spike in this cell
        end
        currentBin = ispk+1;  % Move to next bin
        % --  Update # of samples per iter ---
        muISI = currentBin/(sum(nsp));
        nbinsPerEval = max(20, round(1.5*muISI)); 
    end
end
 response = response(:,maxFilterLength + 1:end);
end