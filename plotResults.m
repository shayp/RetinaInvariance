function plotResults()
load('globalParams.mat');
load('FinalNeuronParameters.mat');
numOfSubPlots = numOfNeurons + 2;
numofRows = ceil(numOfSubPlots / 2);
for i = 1:numOfNeurons
   figure();
    % Plot STA estimator
    subplot(numofRows,2,1);
    plot(NeuronParameters(i).StimulusFilter);
    hold on;
    plot(NeuronParameters(i).expStimulusFilter);
    legend('Learned STA','Experiment STA');
    xlabel('Time before spike');
    ylabel('intensity');
    title(['STA estimator - Neuron: ' num2str(NeuronParameters(i).neuronIndex)]);

    for j = 1:numOfNeurons
        % Plot leaned spike history filter
        subplot(numofRows,2,j + 1);
        plot(NeuronParameters(i).couplingFilters(j,:));
        title(['coupling filter Neuron: ' num2str(NeuronParameters(i).coupledNeurons(j))]);
        xlabel('Time after  spike');
        ylabel('Firing factor');
    end
        lengthOfRepeat = length(NeuronParameters(i).realSpikeRate);
        subplot(numofRows,2,numOfSubPlots);
        plot(16* (1:lengthOfRepeat), NeuronParameters(i).realSpikeRate, 16 * (1:lengthOfRepeat),NeuronParameters(i).simulatedSpikeRate);
        %xlim([100 500]);
        title([ 'Neuron: ' num2str(NeuronParameters(i).neuronIndex) ' R ' num2str(NeuronParameters(i).spikeRateCorrelation)]);
        legend('Experiment firing rate','Simulated firing rate');
        xlabel('Time (ms)');
        ylabel('Firing rate (spikes/sec) ');
end

