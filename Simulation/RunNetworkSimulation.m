function response = RunNetworkSimulation(Stimulus, stimulusFilter, spikeHistoryFilter, meanFiringRate, deltaT, couplingFlag, couplingFilters, coupledSpikesTrain)
spikeHistoryFilterLength = length(spikeHistoryFilter);
if couplingFlag == 1
    numOfCoupledNeurons = size(couplingFilters, 1);
else
    numOfCoupledNeurons = 0;
end
simulationLength = length(Stimulus);
response = zeros(simulationLength, 1);
baseValue = zeros(simulationLength, 1);
historyValue = zeros(simulationLength, 1);

nbinsPerEval = 100;  % Default number of bins to update for each spike

baseValue = baseValue + conv(Stimulus, stimulusFilter, 'same') + meanFiringRate;
for neuronIndex = 1:numOfCoupledNeurons
    coupledValue = conv(coupledSpikesTrain(neuronIndex,:), couplingFilters(neuronIndex,:), 'same')';
    baseValue = baseValue + coupledValue ;
end

rprev = 0;
nsp =0;
tspnext = exprnd(1);
currentBin= 1;
while currentBin <= simulationLength
    iinxt = currentBin:min(currentBin+nbinsPerEval-1,simulationLength);
    rrnxt =  exp(baseValue(iinxt) + historyValue(iinxt)) * deltaT; % Cond Intensity

    rrcum = cumsum(rrnxt) + rprev;  % Cumulative intensity
    if (tspnext >= rrcum(end)) % No spike in this window
            currentBin = iinxt(end)+1;
            rprev = rrcum(end);
    else % Spike!
        ispk =  iinxt(find(rrcum>=tspnext, 1, 'first'));
        nsp = nsp + 1;
        response(ispk) = 1;
        % Record this spike
        mxi = min(simulationLength, ispk+spikeHistoryFilterLength); % determine bins for adding h current
        iiPostSpk = ispk+1:mxi; % time bins affected by post-spike kernel

        if ~isempty(iiPostSpk)
            baseValue(iiPostSpk) = baseValue(iiPostSpk) + spikeHistoryFilter(1:length(iiPostSpk))';
        end
        tspnext = exprnd(1);  % draw next spike time
        rprev = 0; % reset integrated intensity
        currentBin = ispk+1;  % Move to next bin
        % --  Update # of samples per iter ---
        muISI = currentBin/nsp;
        nbinsPerEval = max(20, round(1.5*muISI)); 
    end
end
end