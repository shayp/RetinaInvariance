function plotResults()
load('globalParams.mat');
load('FinalNeuronParameters.mat');
numOfSubPlots = numOfNeurons + 3;
numofRows = ceil(numOfSubPlots / 2);
%'GLM_Full_NeuronParameters', 'GLM_Partial_NeuronParameters', 'LN_NeuronParameters'
for i = 1:numOfNeurons
fig = figure('visible', 'on');
    % Plot STA estimator
    subplot(numofRows,2,1);
    timeBforeSpike = linspace(-0.4,0, 200);
    plot(timeBforeSpike, GLM_Partial_NeuronParameters(i).stimulusFilter,...
        timeBforeSpike, GLM_Full_NeuronParameters(i).stimulusFilter,...
        timeBforeSpike, LN_NeuronParameters(i).stimulusFilter,...
        timeBforeSpike, NeuronParameters(i).expStimulusFilter);

    legend('Partial GLM filter', 'Full GLM filter','LN filter','Experiment filter');
    xlabel('Time before spike(s)');
    ylabel('intensity');
    title(['Stimulus filter estimator - Neuron: ' num2str(NeuronParameters(i).neuronIndex)]);

    for j = 1:numOfNeurons
        % Plot leaned spike history filter
        subplot(numofRows,2,j + 1);
        binsOfCouplingLength = length(GLM_Full_NeuronParameters(i).couplingFilters(j,:));
        timeAfterSpike = linspace(0, binsOfCouplingLength / 500, binsOfCouplingLength);
        plot(timeAfterSpike, GLM_Partial_NeuronParameters(i).couplingFilters(j,:),...
             timeAfterSpike, GLM_Full_NeuronParameters(i).couplingFilters(j,:));
         legend('Partial GLM', 'Full GLM ');
        title(['coupling filter Neuron: ' num2str(NeuronParameters(i).coupledNeurons(j))]);
        xlabel('Time after  spike(s)');
        ylabel('Firing factor');
    end
        lengthOfRepeat = length(NeuronParameters(i).realSpikeRate);
        
        subplot(numofRows,2,numOfSubPlots - 1);
        x = [min(NeuronParameters(i).correaltionVector(1,:)) max(NeuronParameters(i).correaltionVector(1,:))];
        plot(x, x, '-', NeuronParameters(i).correaltionVector(1,:), NeuronParameters(i).correaltionVector(2,:), 'bo',...
            NeuronParameters(i).correaltionVector(1,:), NeuronParameters(i).correaltionVector(3,:), 'ro',...
            NeuronParameters(i).correaltionVector(1,:), NeuronParameters(i).correaltionVector(4,:), 'go');
        legend(' ','Partial GLM filter', 'Full GLM filter','LN filter');
        title([ 'Neuron: ' num2str(NeuronParameters(i).neuronIndex) ' Firing rate correlation ']);
        xlabel('Experiment firing rate');
        ylabel('Simulated firing rate');
        
        subplot(numofRows,2,numOfSubPlots);
        binSize = 1 / 500 * 8;
        spikeRateTime = linspace(0, lengthOfRepeat * binSize, lengthOfRepeat);
        plot(spikeRateTime, NeuronParameters(i).realSpikeRate,...
             spikeRateTime,NeuronParameters(i).glmPartialSimulatedSpikeRate,...
             spikeRateTime,NeuronParameters(i).glmFullSimulatedSpikeRate,...
             spikeRateTime,NeuronParameters(i).lnSimulatedSpikeRate);
        %xlim([100 500]);
        title([ 'Neuron: ' num2str(NeuronParameters(i).neuronIndex)]);
        legend('Experiment firing rate',['GLM partial R = ' num2str(NeuronParameters(i).spikeRateCorrelation(1))],...
            ['GLM Full R = ' num2str(NeuronParameters(i).spikeRateCorrelation(2))],...
            ['LN R = ' num2str(NeuronParameters(i).spikeRateCorrelation(3))]);
        xlabel('Time (s)');
        ylabel('Firing rate (spikes/sec) ');
        savefig(fig,['./Graphs/Neuron_' num2str(NeuronParameters(i).neuronIndex) '_Results']);
end

