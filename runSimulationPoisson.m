function response = runSimulationPoisson(numOfNeurons, Stimulus, Filters, stimulusFilterLength, couplingFilterLength,deltaT, fCoupling)

maxFilterLength = max(stimulusFilterLength, couplingFilterLength);
simulationLength = length(Stimulus);
response = zeros(numOfNeurons, simulationLength);
projected = zeros(numOfNeurons, simulationLength);
nCount = 0;
for neuronIndex = 1:numOfNeurons
    projected(neuronIndex, :) = conv(Stimulus, Filters(neuronIndex).stimulusFilter, 'same') + Filters(neuronIndex).meanFiringRate;
end
for i = 1:simulationLength - 1
        curentLambdas = exp(projected(:,i)) * deltaT;
        samples = poissrnd(curentLambdas);
        spikedNeurons = find(samples);
        response(spikedNeurons, i) = 1;
        if fCoupling == 1 && isempty(spikedNeurons) == 0
            minStep = min(couplingFilterLength,simulationLength - i);
            projected(1:numOfNeurons, i + 1: i + minStep) = projected(1:numOfNeurons, i + 1:i + minStep) + Filters(1:numOfNeurons).couplingFilters(spikedNeurons,1:minStep);
        end
end
end
