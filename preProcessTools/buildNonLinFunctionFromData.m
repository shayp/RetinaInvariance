function [sigCurveParams, xAxis,yAxis] = buildNonLinFunctionFromData(binnedStimulusConv, binnedSpikes, numofParameters)
length(binnedSpikes)
length(binnedStimulusConv)
lengthOfExp = min(length(binnedSpikes), length(binnedStimulusConv));
maxStimulusConv = max(binnedStimulusConv);
minStimulusConv = min(binnedStimulusConv);

binsValues = linspace(minStimulusConv, maxStimulusConv, numofParameters);
hitNumner = zeros(numofParameters, 1);
binnedStimulusSum = zeros(numofParameters, 1);
binnedSpikesSum = zeros(numofParameters, 1);

for i = 1:lengthOfExp
    [~, wantedIndex] = min(abs(binsValues - binnedStimulusConv(i)));
    hitNumner(wantedIndex) = hitNumner(wantedIndex) + 1;
    binnedStimulusSum(wantedIndex) = binnedStimulusSum(wantedIndex) + binnedStimulusConv(i);
    binnedSpikesSum(wantedIndex) = binnedSpikesSum(wantedIndex) + binnedSpikes(i);
end

for i = 1:numofParameters
    if hitNumner(i) ~= 0
        binnedStimulusSum(i) = binnedStimulusSum(i) / hitNumner(i);
        binnedSpikesSum(i) = binnedSpikesSum(i) / hitNumner(i);
    end
end

cumlativeOfPoints = cumsum(hitNumner);
maxIndexes = find(cumlativeOfPoints > 0.9 * lengthOfExp);

initValue = [0, 0.1, 1, 1, 1, 0];
[xAxis,yAxis,sigCurveParams] = fitSigmoidToData(binnedStimulusSum, binnedSpikesSum, -2,binnedStimulusSum(maxIndexes(1)), initValue);

% figure();
% plot(binnedStimulusSum, binnedSpikesSum,'.', binnedStimulusSum, sigCurveParams(binnedStimulusSum));
% title('Non linear function');
% legend('Data', 'Fit');
end