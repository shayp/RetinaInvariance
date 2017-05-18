load SpikeTimeGaussian

load GaussianMovieData
%%

K=find(WW(:,1)>150);StimWW=K(find(diff(K)>10));

Stim=zeros(601000*20,1);
NewStimTime=[];

for i=106:124% 1:20 %length(StimWW)
    i
    NewStimTime=[NewStimTime StimTime(i)-StimTime(106)];
    for j=3:1790
        % Stim(round((j-1)*333.333+(1:334)),1)=WW(StimWW(i)+j,1);
        Stim(StimTime(i)-StimTime(106)+1+round((j-1)*333.333+(1:334)))=WW(StimWW(i)+j,3);        
    end
end


for N=1:100
    tmp=TT(N).sp;
    tmp=tmp(find(tmp>StimTime(106)))-StimTime(106);
    TTT(N).sp=tmp;    
end
%%
STAALL=zeros(5000,10);
figure;hold on
for N=1:50
    N
    STA=zeros(5000,1);
    tmp=TTT(N).sp;
    tmp=tmp(find(tmp>4100 & tmp<(601000*20-4000)));
    for i=1:length(tmp);
        STA=STA+Stim(tmp(i)+(-4000:999));
    end    
    STAALL(:,N)=STA;
    STA=STA-mean(STA);STA=STA/norm(STA);
    plot(STA)
    title(N)
    pause(.1)
end

  
%%

% TTDataForStudents=TTT([1 2 3 26]);
% save CourseData TTDataForStudents Stim NewStimTime
%%

clear 
load CourseData

STAALL=zeros(5000,10);
figure;hold on
for N=1:length(TTDataForStudents)
    N
    STA=zeros(5000,1);
    tmp=TTDataForStudents(N).sp;
    tmp=tmp(find(tmp>4100 & tmp<(601000*20-4000)));
    for i=1:length(tmp);
        STA=STA+Stim(tmp(i)+(-4000:999));
    end    
    STAALL(:,N)=STA;
    STA=STA-mean(STA);STA=STA/norm(STA);
    plot(STA)
    title(N)
    pause(.1)
end
