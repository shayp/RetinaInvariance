function response = RunSimulation(numOfNeurons, Stimulus, neuronParameters, stimulusFilterLength, deltaT, couplingFlag)
if couplingFlag == 1
    couplingFilterLength = size(neuronParameters(1).couplingFilters,2);
else
    couplingFilterLength = 1;
end
simulationLength = length(Stimulus);
response = zeros(numOfNeurons, simulationLength);
baseValue = zeros(numOfNeurons, simulationLength);
couplingValue = zeros(numOfNeurons, simulationLength);

nCount = 0;
nbinsPerEval = 100;  % Default number of bins to update for each spike

for neuronIndex = 1:numOfNeurons
    baseValue(neuronIndex,:) = conv(Stimulus, neuronParameters(neuronIndex).stimulusFilter, 'same');
    baseValue(neuronIndex,:) = baseValue(neuronIndex,:) + neuronParameters(neuronIndex).meanFiringRate ;
end

rprev = zeros(1,numOfNeurons);
nsp = zeros(1,numOfNeurons);
tspnext = exprnd(1,1,numOfNeurons);
currentBin= 1;
while currentBin <= simulationLength
    iinxt = currentBin:min(currentBin+nbinsPerEval-1,simulationLength);
    nii = length(iinxt);  % Number of bins
    rrnxt =  exp(baseValue(:,iinxt) + couplingValue(:,iinxt)) * deltaT; % Cond Intensity

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
        mxi = min(simulationLength, ispk+couplingFilterLength); % determine bins for adding h current
        iiPostSpk = ispk+1:mxi;
         for ic = 1:length(spcells)
            icell = spcells(ic);
            nsp(icell) = nsp(icell)+1;
            response(icell, ispk) = 1;
            %couplingValue(icell,iiPostSpk) = zeros(1,length(iiPostSpk));
            if couplingFlag == 1
                for i = 1:numOfNeurons
                    couplingValue(i,iiPostSpk) = couplingValue(i,iiPostSpk) + neuronParameters(i).couplingFilters(icell,1:length(iiPostSpk));
                end
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
end