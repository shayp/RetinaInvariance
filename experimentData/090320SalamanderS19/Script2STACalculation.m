load SpikeTimeFatTailed

load FatTailedMovieData

%%

STAAll=zeros(5000,10);

for N=1:10
    STA=zeros(5000,1);
    K=find(WW(:,1)>150);StimWW=K(find(diff(K)>10));

    for i=1:20 %length(StimWW)
        [N i]
        Stim=zeros(610000,2);
        for j=3:1790
            Stim(round((j-1)*333.333+(1:334)),1)=WW(StimWW(i)+j,1);
            Stim(round((j-1)*333.333+(1:334)),2)=WW(StimWW(i)+j,3);        
        end
    %     plot(E(StimTime(i):(StimTime(i)+20000)));hold on;
    %     plot(Stim(1:20000)*10+500,'.-r')
    %     pause
    %     hold off
        tmp=TT(N).sp;
        tmp=tmp(find(tmp>(StimTime(i)+4500) & tmp<StimTime(i+1)))-StimTime(i)+1;
        for j=1:length(tmp)
            STA=STA+Stim((tmp(j)-4000):(tmp(j)+999),2);
        end
    end
    
    STAAll(:,N)=STA;
end

    
%%    
figure;hold on
for i=1:10
    tmp=STAAll(:,i)-mean(STAAll(:,i));
    tmp=tmp/max(abs(tmp));
    plot(((1:5000)-4000)/10,tmp)    
end

