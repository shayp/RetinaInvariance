function response = runSimulationPoisson(numOfNeurons, Stimulus, Filters, stimulusFilterLength, couplingFilterLength,deltaT)

maxFilterLength = max(stimulusFilterLength, couplingFilterLength);
simulationLength = length(Stimulus);
response = zeros(numOfNeurons, simulationLength + maxFilterLength);
projected = zeros(numOfNeurons, simulationLength + maxFilterLength);
nCount = 0;
for neuronIndex = 1:numOfNeurons
    projected(neuronIndex, stimulusFilterLength + 1:end) = conv(Stimulus, Filters(neuronIndex).StimulusFilter, 'same') + Filters(neuronIndex).meanFiringRate;
end
for i = maxFilterLength + 1:simulationLength
    for neuronIndex = 1:numOfNeurons
        curentLambda = exp(projected(neuronIndex,i)) * deltaT;
        sample = poissrnd(curentLambda);
        if sample >= 1
            nCount  = nCount + 1;
            %sample = 1;
            for couplingIndex = 1:numOfNeurons
                projected(couplingIndex, i + 1: i + couplingFilterLength) = projected(couplingIndex, i + 1: i + couplingFilterLength) + Filters(couplingIndex).couplingFilters(neuronIndex,:);
            end
        end
      response(neuronIndex, i) = sample;
        
    end
end
response = response(:,maxFilterLength + 1:end);

end