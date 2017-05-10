% Load the data for GLM
datdir = './';  
load([datdir, 'Stim']);    
load([datdir,'stimtimes']); 
load([datdir, 'SpTimes']); 
ncells = length(SpTimes);
couplenNeurons = [2];
runGLM(1, Stim, stimtimes, SpTimes, couplenNeurons)
