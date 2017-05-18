 ihbasprs.ncols = 5; 
 
 ihbasprs.hpeaks = [.1 2];  
 ihbasprs.b = .5;  
 ihbasprs.absref = .1;  %% (optional)
 dt = 0.01;
 [iht,ihbas,ihbasis] = CosineBumps(ihbasprs,dt);
 
 figure();
 subplot(2,1,1);
 plot(ihbasis);
 subplot(2,1,2);
 plot(ihbas(:,3) -1 * ihbas(:,1));