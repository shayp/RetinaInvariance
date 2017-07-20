
clear all;
datdir = './';  
load([datdir, 'Stim']);    
load([datdir,'stimtimes']); 
load([datdir, 'SpTimes']); 

lastIndex = stimtimes(end);
wantedSampFactor = 20;
tmp = changeSpikeRsolution(SpTimes, lastIndex, wantedSampFactor)