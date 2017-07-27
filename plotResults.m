function plotResults(numOfNeurons)
load('globalParams.mat');
load('FinalNeuronParameters.mat');
numOfSubPlots = numOfNeurons + 5;
numofRows = ceil(numOfSubPlots / 2);
for i = 1:numOfNeurons
fig = figure('visible', 'on');

    % Plot stimulus filter
    subplot(numofRows,2,1);
    timeBforeSpike = linspace(- stimulusFilterSizeForSimulation / binsInSecond,0, stimulusFilterSizeForSimulation)';

    plot(timeBforeSpike, NeuronParameters(i).stimulusFilter,...
         timeBforeSpike, lnParams(i).stimulusFilter,...
         timeBforeSpike, glmPartialParams(i).stimulusFilter,...
         timeBforeSpike, glmFullParams(i).stimulusFilter);

    legend('STA','LN ', 'Partial GLM ', 'Full GLM ','location', 'bestoutside');
    xlabel('Time before spike (s)');
    ylabel('intensity');
    title(['Stimulus filter estimator - Neuron: ' num2str(NeuronParameters(i).neuronIndex)]);
    
    X = -10:0.1:2;
    Y = exp(X) * deltaT;
    % Plot non lineartity
    subplot(numofRows,2,2);
    plot(X, exp(lnParams(i).meanFiringRate) * Y,...
         X, exp(glmPartialParams(i).meanFiringRate) * Y,...
         X, exp(glmFullParams(i).meanFiringRate) * Y,...
         X, lnBusgangParams(i).expFunction(X));

    legend('LN','Partial GLM', 'Full GLM', 'LN Busgang','location', 'bestoutside');
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
        xlabel('Time after  spike (s)');
        ylabel('Firing factor');
    end
        lengthOfRepeat = length(NeuronParameters(i).realSpikeRate);
        
        subplot(numofRows,2,numOfSubPlots - 2);
        x = [min(NeuronParameters(i).correaltionVector(1,:)) max(NeuronParameters(i).correaltionVector(1,:))];
        plot(x, x, '-',...
             NeuronParameters(i).correaltionVector(1,:), NeuronParameters(i).correaltionVector(2,:), 'bo',...
             NeuronParameters(i).correaltionVector(1,:), NeuronParameters(i).correaltionVector(3,:), 'ro',...
             NeuronParameters(i).correaltionVector(1,:), NeuronParameters(i).correaltionVector(4,:), 'go',...
             NeuronParameters(i).correaltionVector(1,:), NeuronParameters(i).correaltionVector(5,:), 'ko',...
             NeuronParameters(i).correaltionVector(1,:), NeuronParameters(i).correaltionVector(6,:), 'mo');

        legend('y = x','LN ','Partial GLM ', 'Full GLM', 'LN Busgang', '1/2 of experiment','location', 'bestoutside');
        title('Firing rate correlation');
        xlabel('Experiment firing rate');
        ylabel('Simulated firing rate');

        subplot(numofRows,2,numOfSubPlots - 1);
        binSize = deltaT * windowSizeForFiringRate;
        spikeRateTime = linspace(0, lengthOfRepeat * binSize, lengthOfRepeat);
        plot(spikeRateTime, NeuronParameters(i).realSpikeRate,...
             spikeRateTime,NeuronParameters(i).lnSimulatedSpikeRate,...
             spikeRateTime,NeuronParameters(i).glmPartialSimulatedSpikeRate,...
             spikeRateTime,NeuronParameters(i).glmFullSimulatedSpikeRate,...
             spikeRateTime,NeuronParameters(i).lnBusgangSpikeRate);

         title('Firing rate');
        legend('Experiment firing rate',...
            ['LN R = ' num2str(NeuronParameters(i).spikeRateCorrelation(1))],...
            ['GLM partial R = ' num2str(NeuronParameters(i).spikeRateCorrelation(2))],...
            ['GLM Full R = ' num2str(NeuronParameters(i).spikeRateCorrelation(3))],...
            ['LN Busgang R = ' num2str(NeuronParameters(i).spikeRateCorrelation(4))],...
            'location', 'bestoutside');
        xlabel('Time (s)');
        ylabel('Firing rate (spikes/sec) ');
        
        subplot(numofRows,2,numOfSubPlots);
        strNames = {'LN','GLM Partial','GLM Full', 'LN Busgang'};
        bar(NeuronParameters(i).perecentExplained)
        title('% Explained from 1/2 of experiment');
        set(gca, 'XTickLabel', strNames, 'XTick', 1:numel(strNames));
        
        savefig(fig,['./Graphs/Neuron_' num2str(NeuronParameters(i).neuronIndex) '_Results']);

end

