function response = RunGLMSimulation(numOfNeurons, Stimulus, Filters, stimulusFilterLength, couplingFilterLength,deltaT)

maxFilterLength = max(stimulusFilterLength, couplingFilterLength);
simulationLength = length(Stimulus);
response = zeros(numOfNeurons, simulationLength + maxFilterLength);
projectedStimulus = zeros(numOfNeurons, simulationLength + maxFilterLength);
nCount = 0;
for neuronIndex = 1:numOfNeurons
    projectedStimulus(neuronIndex, stimulusFilterLength + 1:end) = conv(Stimulus, Filters(neuronIndex).StimulusFilter, 'same');
end
for i = maxFilterLength + 1:simulationLength
    for neuronIndex = 1:numOfNeurons

        projectionTrm = projectedStimulus(neuronIndex,i)+ Filters(neuronIndex).meanFiringRate;
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