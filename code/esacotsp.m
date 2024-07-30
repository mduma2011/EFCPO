
function [res] = esacotsp(paraminfo)

  % Effective Strategies Ant Colony Optimization (ESACO). This algorithm was written as part of the
  % paper "Effective Heuristics for Ant Colony Optimization to Handle Large-Scale Problems"
  % published in 2017. We need to see if we can improve the performance of
  % the algorithm using various inputs
  % By default, the localsearch is 2-opt. for now.
  
   antpop = paraminfo.antpop;
   maxiteration = paraminfo.maxiteration;   
   alpha = paraminfo.alpha;
   rho = paraminfo.rho;
   beta = paraminfo.beta;  
   maxclcount = paraminfo.maxclcount;
   numberofcandidatesconst = paraminfo.numberofcandidatesconst;
   numberofcandidates =  paraminfo.numberofcandidates;
   fullpathname  =  paraminfo.fullpathname; 
   
   [BESTANT0, ETA,bestant, besttour] = esacomain(antpop,maxiteration,alpha,rho,beta,maxclcount,numberofcandidatesconst,numberofcandidates, fullpathname);
   
   res.bestant = bestant;
   res.besttour = besttour;
   res.ETA = ETA;  
   res.ants = BESTANT0;   
end

