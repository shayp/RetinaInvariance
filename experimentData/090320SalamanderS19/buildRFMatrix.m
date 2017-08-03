function buildRFMatrix()
load ('SpikeTimeGaussian');
numOfNeurons = length(TT);
TT
figure();
for i = 1:numOfNeurons
    RFMatrix(i).data = AmpAndPinToLocation(TT(i).Amp,TT(i).channel);
    hold on;
    plot(RFMatrix(i).data(1), RFMatrix(i).data(2), '.');
    text(RFMatrix(i).data(1), RFMatrix(i).data(2) , num2str(i));
    hold off;
    drawnow;
end

save('RFMatrix','RFMatrix');
end