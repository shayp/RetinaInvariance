function bla(stimulus,spikes, lengthOfSTA)
STA = zeros(lengthOfSTA,1);
spikes(1:5)
numOfSPikes = length(spikes)
for i = 1:numOfSPikes
    spikes(i)
    STA = STA + stimulus(spikes(i) - lengthOfSTA:spikes(i) - 1);
end
STA = STA / numOfSPikes;
figure();
plot(STA);
end