numOfNeruons = 251;
neuronCorretlation = zeros(numOfNeruons, 1);
load('RFMatrix');
pathFolder = './Results/Results_1_8_17-HistoryAll/';
d = dir(pathFolder);
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

numOfFolders = length(nameFolds);

for i = 1:numOfFolders
load([pathFolder nameFolds{i} '/FinalNeuronParameters']);
    neuronCorretlation(NeuronParameters(1).neuronIndex)  = NeuronParameters(1).spikeRateCorrelation(1);
end
figure();
for i = 1:numOfNeruons
hold on;
    if neuronCorretlation(i) > 0.6
        plot(RFMatrix(i).data(1), RFMatrix(i).data(2), '.r');
        text(RFMatrix(i).data(1), RFMatrix(i).data(2) , [num2str(i) ', ' num2str(neuronCorretlation(i),2)]);
    elseif neuronCorretlation(i) > 0.5
       plot(RFMatrix(i).data(1), RFMatrix(i).data(2), '.g');
       text(RFMatrix(i).data(1), RFMatrix(i).data(2) , [num2str(i) ', ' num2str(neuronCorretlation(i),2)]);
    elseif neuronCorretlation(i) > 0.4
       plot(RFMatrix(i).data(1), RFMatrix(i).data(2), '.b');
       text(RFMatrix(i).data(1), RFMatrix(i).data(2) , [num2str(i) ', ' num2str(neuronCorretlation(i),2)]);
     else
       plot(RFMatrix(i).data(1), RFMatrix(i).data(2), '.k');
     end

    hold off;
    drawnow;
end

