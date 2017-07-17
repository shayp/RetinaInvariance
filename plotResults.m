function plotResults()
load('globalParams.mat');
load('FinalNeuronParameters.mat');
numOfSubPlots = numOfNeurons + 3;
numofRows = ceil(numOfSubPlots / 2);
for i = 1:numOfNeurons
fig = figure('visible', 'off');
    % Plot STA estimator
    subplot(numofRows,2,1);
    timeBforeSpike = linspace(-0.4,0, 200);
    plot(timeBforeSpike, NeuronParameters(i).StimulusFilter);
    hold on;
    plot(timeBforeSpike, NeuronParameters(i).expStimulusFilter);
    legend('Learned stimulus filter','Experiment stimulus filter');
    xlabel('Time before spike(s)');
    ylabel('intensity');
    title(['Stimulus filter estimator - Neuron: ' num2str(NeuronParameters(i).neuronIndex)]);

    for j = 1:numOfNeurons
        % Plot leaned spike history filter
        subplot(numofRows,2,j + 1);
        binsOfCouplingLength = length(NeuronParameters(i).couplingFilters(j,:));
        timeAfterSpike = linspace(0, binsOfCouplingLength / 500, binsOfCouplingLength);
        plot(timeAfterSpike, NeuronParameters(i).couplingFilters(j,:));
        title(['coupling filter Neuron: ' num2str(NeuronParameters(i).coupledNeurons(j))]);
        xlabel('Time after  spike(s)');
        ylabel('Firing factor');
    end
        lengthOfRepeat = length(NeuronParameters(i).realSpikeRate);
        
        subplot(numofRows,2,numOfSubPlots - 1);
        x = [min(NeuronParameters(i).correaltionVector(1,:)) max(NeuronParameters(i).correaltionVector(1,:))];
        plot(x, x, '-', NeuronParameters(i).correaltionVector(1,:), NeuronParameters(i).correaltionVector(2,:), 'bo');
        title([ 'Neuron: ' num2str(NeuronParameters(i).neuronIndex) ' Firing rate correlation ' num2str(NeuronParameters(i).spikeRateCorrelation)]);
        xlabel('Experiment firing rate');
        ylabel('Simulated firing rate');
        
        subplot(numofRows,2,numOfSubPlots);
        plot(16* (1:lengthOfRepeat), NeuronParameters(i).realSpikeRate, 16 * (1:lengthOfRepeat),NeuronParameters(i).simulatedSpikeRate);
        %xlim([100 500]);
        title([ 'Neuron: ' num2str(NeuronParameters(i).neuronIndex) ' R ' num2str(NeuronParameters(i).spikeRateCorrelation)]);
        legend('Experiment firing rate','Simulated firing rate');
        xlabel('Time (ms)');
        ylabel('Firing rate (spikes/sec) ');
        savefig(fig,['./Graphs/Neuron_' num2str(NeuronParameters(i).neuronIndex) '_Results']);
end

