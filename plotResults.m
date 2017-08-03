function plotResults(numOfNeurons, filePath)
load('globalParams.mat');
load([filePath 'FinalNeuronParameters']);
numOfSubPlotsFig1 = numOfNeurons + 2;
numofRowsFig1 = ceil(numOfSubPlotsFig1 / 2);
numOfSubPlotsFig2 = 4;
numofRowsFig2 = ceil(numOfSubPlotsFig2 / 2);
mumOfRepeatsToPlot = 10;
for i = 1:numOfNeurons
    fig1 = figure('visible', 'on');

    % Plot stimulus filter
    subplot(numofRowsFig1,2,1);
    timeBforeSpike = linspace(- stimulusFilterSizeForSimulation / binsInSecond,0, stimulusFilterSizeForSimulation)';

    plot(timeBforeSpike, NeuronParameters(i).stimulusFilter,...
         timeBforeSpike, lnParams(i).stimulusFilter,...
         timeBforeSpike, glmFullParams(i).stimulusFilter,...
         timeBforeSpike, NeuronParameters(i).repSTA);

    legend('Non repeat STA','LN ', 'Full GLM ','Repeat STA','location', 'bestoutside');
    xlabel('Time before spike (s)');
    ylabel('intensity');
    title(['Stimulus filter estimator - Neuron: ' num2str(NeuronParameters(i).neuronIndex)]);
    
    X = -10:0.1:2.5;
    Y = exp(X) * deltaT;
    % Plot non lineartity
    subplot(numofRowsFig1,2,2);
    plot(X, exp(lnParams(i).meanFiringRate) * Y,...
         X, exp(glmFullParams(i).meanFiringRate) * Y,...
         X, lnBusgangParams(i).expFunction(X));

    legend('LN', 'Full GLM', 'LN Busgang','location', 'bestoutside');
    xlabel('X');
    ylabel('exp(X)');
    title(['nonLinear Exp: ' num2str(NeuronParameters(i).neuronIndex)]);
    
    for j = 1:numOfNeurons
        % Plot leaned spike history filter
        subplot(numofRowsFig1,2, j + 2);
        binsOfCouplingLength = length(glmFullParams(i).couplingFilters(j,:));
        timeAfterSpike = linspace(0, binsOfCouplingLength / binsInSecond, binsOfCouplingLength);
        plot(timeAfterSpike, exp(glmFullParams(i).couplingFilters(j,:)),...
             timeAfterSpike, ones(1, length(timeAfterSpike)), '--');
         legend('Full GLM','location', 'bestoutside');
        title(['coupling filter Neuron: ' num2str(NeuronParameters(i).neuronsInNetwork(j))]);
        xlabel('Time after  spike (s)');
        ylabel('Firing factor');
    end
    
        savefig(fig1,[filePath 'Neuron_' num2str(NeuronParameters(i).neuronIndex) '_Filters']);
        fig2 = figure('visible', 'on');

        lengthOfRepeat = length(NeuronParameters(i).realSpikeRate);
        
        subplot(numofRowsFig2,2,1);
        x = [min(NeuronParameters(i).correaltionVector(1,:)) max(NeuronParameters(i).correaltionVector(1,:))];
        plot(x, x, '-',...
             NeuronParameters(i).correaltionVector(1,:), NeuronParameters(i).correaltionVector(2,:), 'bo',...
             NeuronParameters(i).correaltionVector(1,:), NeuronParameters(i).correaltionVector(3,:), 'go',...
             NeuronParameters(i).correaltionVector(1,:), NeuronParameters(i).correaltionVector(4,:), 'ko',...
             NeuronParameters(i).correaltionVector(1,:), NeuronParameters(i).correaltionVector(5,:), 'mo');

        legend('y = x','LN ', 'Full GLM', 'LN Busgang', '1/2 of experiment','location', 'bestoutside');
        title('Firing rate correlation');
        xlabel('Experiment firing rate');
        ylabel('Simulated firing rate');

        subplot(numofRowsFig2, 2, 2);
        binSize = deltaT * windowSizeForFiringRate;
        spikeRateTime = linspace(0, lengthOfRepeat * binSize, lengthOfRepeat);
        startLimit = randi(lengthOfRepeat - 201);
        plot(spikeRateTime, NeuronParameters(i).realSpikeRate,...
             spikeRateTime,NeuronParameters(i).lnSimulatedSpikeRate,...
             spikeRateTime,NeuronParameters(i).glmFullSimulatedSpikeRate,...
             spikeRateTime,NeuronParameters(i).lnBusgangSpikeRate);

        title('Firing rate');
        legend(['Experiment firing rate. R = ' num2str(NeuronParameters(i).spikeRateCorrelation(4))],...
            ['LN R = ' num2str(NeuronParameters(i).spikeRateCorrelation(1))],...
            ['GLM Full R = ' num2str(NeuronParameters(i).spikeRateCorrelation(2))],...
            ['LN Busgang R = ' num2str(NeuronParameters(i).spikeRateCorrelation(3))],...
            'location', 'bestoutside');
        xlabel('Time (s)');
        ylabel('Firing rate (spikes/sec) ');
        xlim([spikeRateTime(startLimit) spikeRateTime(startLimit + 200)]);
        
%         subplot(numofRowsFig2,2,3);
%         strNames = {'LN','GLM Full', 'LN Busgang'};
%         bar(NeuronParameters(i).perecentExplained)
%         title(['% R Explained from 1/2 of experiment. R = ' num2str(NeuronParameters(i).spikeRateCorrelation(4))]);
%         set(gca, 'XTickLabel', strNames, 'XTick', 1:numel(strNames));
        
         subplot(numofRowsFig2,2,3);
         for j = 1:mumOfRepeatsToPlot
             realRaster = find(NeuronParameters(i).scaledRepSpikes(j,:));
             ExpPoints = zeros(1, length(realRaster));
             ExpPoints = ExpPoints + 20 + j;
             hold on;
             plot (realRaster * deltaT, ExpPoints, '.k');
             hold off;
             
             LNRaster = find(NeuronParameters(i).LNSimulation(j,:));
             LNPoints = zeros(1, length(LNRaster));
             LNPoints = LNPoints + 10 + j;
             hold on;
             plot (LNRaster * deltaT, LNPoints, '.r');
             hold off;
             
             GLMFullRaster = find(NeuronParameters(i).GLMFullSimulation(j,:));
             GLMFullPoints = zeros(1, length(GLMFullRaster));
             GLMFullPoints = GLMFullPoints  + j;
             hold on;
             plot (GLMFullRaster * deltaT, GLMFullPoints, '.b');
             hold off;
         end
         %xlim([(startLimit * binSize) (startLimit * binSize + 200 * binSize)]);
         xlim([spikeRateTime(startLimit) spikeRateTime(startLimit + 200)]);
         set(gca, 'YTickLabel', '');
         xlabel('Time (s)');
         title('Raster Plot');
         legend('Experiment data', 'LN Simulation','Full GLM Simulation','location', 'bestoutside');

        subplot(numofRowsFig2,2,4);
        [histRep, edgesRep] = histcounts(NeuronParameters(i).repISI);
        histRep = histRep / length(NeuronParameters(i).repISI);
        histRep = [histRep histRep(end)];
        hold on;
        stairs(edgesRep, histRep);
        hold off;
        
        [histLN, edgesLN] = histcounts(NeuronParameters(i).LNISI);
        histLN = histLN / length(NeuronParameters(i).LNISI);
        histLN = [histLN histLN(end)];
        hold on;
        stairs(edgesLN, histLN);
        hold off;
        
        
        [histGLMFull, edgesGLMFull] = histcounts(NeuronParameters(i).GLMFullISI);
        histGLMFull = histGLMFull / length(NeuronParameters(i).GLMFullISI);
        histGLMFull = [histGLMFull histGLMFull(end)];
        hold on;
        stairs(edgesGLMFull, histGLMFull);
        hold off;
         
         legend('Experiment','LN ', 'Full GLM','location', 'bestoutside');
         title('Inter spike interval');
        savefig(fig2,[filePath 'Neuron_' num2str(NeuronParameters(i).neuronIndex) '_Results']);
end

