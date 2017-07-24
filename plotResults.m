function plotResults(numOfNeurons)
load('globalParams.mat');
load('FinalNeuronParameters.mat');
numOfSubPlots = numOfNeurons + 4;
numofRows = ceil(numOfSubPlots / 2);
%'glmFullParams', 'glmPartialParams', 'lnParams'
for i = 1:numOfNeurons
fig = figure('visible', 'on');

    % Plot stimulus filter
    subplot(numofRows,2,1);
    timeBforeSpike = linspace(- stimulusFilterSizeForSimulation / binsInSecond,0, stimulusFilterSizeForSimulation)';

    plot(timeBforeSpike, glmPartialParams(i).stimulusFilter,...
        timeBforeSpike, glmFullParams(i).stimulusFilter,...
        timeBforeSpike, lnParams(i).stimulusFilter,...
        timeBforeSpike, NeuronParameters(i).stimulusFilter);

    legend('Partial GLM filter', 'Full GLM filter','LN filter','Experiment filter','location', 'bestoutside');
    xlabel('Time before spike(s)');
    ylabel('intensity');
    title(['Stimulus filter estimator - Neuron: ' num2str(NeuronParameters(i).neuronIndex)]);
    
    X = -2:0.1:1;
    Y = exp(X);
    % Plot non lineartity
    subplot(numofRows,2,2);
    plot(X, exp(glmPartialParams(i).meanFiringRate) * Y,...
        X, exp(glmFullParams(i).meanFiringRate) * Y,...
        X, exp(lnParams(i).meanFiringRate) * Y);

    legend('Partial GLM', 'Full GLM','LN','location', 'bestoutside');
    xlabel('X');
    ylabel('exp(X)');
    title(['nonLinear Exp: ' num2str(NeuronParameters(i).neuronIndex)]);
    
    for j = 1:numOfNeurons
        % Plot leaned spike history filter
        subplot(numofRows,2, j + 2);
        binsOfCouplingLength = length(glmFullParams(i).couplingFilters(j,:));
        timeAfterSpike = linspace(0, binsOfCouplingLength / binsInSecond, binsOfCouplingLength);
        plot(timeAfterSpike, glmPartialParams(i).couplingFilters(j,:),...
             timeAfterSpike, glmFullParams(i).couplingFilters(j,:));
         legend('Partial GLM', 'Full GLM','location', 'bestoutside');
        title(['coupling filter Neuron: ' num2str(NeuronParameters(i).neuronsInNetwork(j))]);
        xlabel('Time after  spike(s)');
        ylabel('Firing factor');
    end
        lengthOfRepeat = length(NeuronParameters(i).realSpikeRate);
        
        subplot(numofRows,2,numOfSubPlots - 1);
        x = [min(NeuronParameters(i).correaltionVector(1,:)) max(NeuronParameters(i).correaltionVector(1,:))];
        plot(x, x, '-', NeuronParameters(i).correaltionVector(1,:), NeuronParameters(i).correaltionVector(2,:), 'bo',...
            NeuronParameters(i).correaltionVector(1,:), NeuronParameters(i).correaltionVector(3,:), 'ro',...
            NeuronParameters(i).correaltionVector(1,:), NeuronParameters(i).correaltionVector(4,:), 'go');
        legend(' ','Partial GLM filter', 'Full GLM filter','LN filter','location', 'bestoutside');
        title([ 'Neuron: ' num2str(NeuronParameters(i).neuronIndex) ' Firing rate correlation ']);
        xlabel('Experiment firing rate');
        ylabel('Simulated firing rate');

        subplot(numofRows,2,numOfSubPlots);
        binSize = 1 / 1000 * windowSizeForFiringRate;
        spikeRateTime = linspace(0, lengthOfRepeat * binSize, lengthOfRepeat);
        plot(spikeRateTime, NeuronParameters(i).realSpikeRate,...
             spikeRateTime,NeuronParameters(i).glmPartialSimulatedSpikeRate,...
             spikeRateTime,NeuronParameters(i).glmFullSimulatedSpikeRate,...
             spikeRateTime,NeuronParameters(i).lnSimulatedSpikeRate);

         title([ 'Neuron: ' num2str(NeuronParameters(i).neuronIndex)]);
        legend('Experiment firing rate',['GLM partial R = ' num2str(NeuronParameters(i).spikeRateCorrelation(1))],...
            ['GLM Full R = ' num2str(NeuronParameters(i).spikeRateCorrelation(2))],...
            ['LN R = ' num2str(NeuronParameters(i).spikeRateCorrelation(3))],'location', 'bestoutside');
        xlabel('Time (s)');
        ylabel('Firing rate (spikes/sec) ');
        savefig(fig,['./Graphs/Neuron_' num2str(NeuronParameters(i).neuronIndex) '_Results']);
end

