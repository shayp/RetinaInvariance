function response = runSimulationPoisson(numOfNeurons, Stimulus, Filters, stimulusFilterLength, couplingFilterLength,deltaT)

maxFilterLength = max(stimulusFilterLength, couplingFilterLength);
simulationLength = length(Stimulus);
response = zeros(numOfNeurons, simulationLength);
projected = zeros(numOfNeurons, simulationLength);
nCount = 0;
for neuronIndex = 1:numOfNeurons
    projected(neuronIndex, :) = conv(Stimulus, Filters(neuronIndex).StimulusFilter, 'same') + Filters(neuronIndex).meanFiringRate;
end
for i = 1:simulationLength - 1
    for neuronIndex = 1:numOfNeurons
        curentLambda = exp(projected(neuronIndex,i)) * deltaT;
        sample = poissrnd(curentLambda);
        if sample >= 1
            nCount  = nCount + 1;
            %sample = 1;
            for couplingIndex = 1:numOfNeurons
                minStep = min(couplingFilterLength,simulationLength - i);
                projected(couplingIndex, i + 1: i + minStep) = projected(couplingIndex, i + 1:i + minStep) + Filters(couplingIndex).couplingFilters(neuronIndex,1:minStep);
            end
        end
      response(neuronIndex, i) = sample;
        
    end
end
end