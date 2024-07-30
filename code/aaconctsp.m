function [res] = aaconctsp(paraminfo)
  % Adaptive Ant Colony Optimization with Node Clustering (AACO-NC). This algorithm was written as part of the
  % paper "Adaptive Ant Colony Optimization with node clustering applied to the Travelling Salesman Problem"
  % published in 2022. We need to see if we can improve the performance of the algorithm using various inputs.
    
   alpha = paraminfo.alpha;
   beta = paraminfo.beta;
   nnsize = paraminfo.nnsize;
   nE = paraminfo.nE;
   Distance = paraminfo.Distance;
   nearneighbours = paraminfo.nearneighbours;
   maxiteration = paraminfo.maxiteration;
   antpop = paraminfo.antpop;
   Q = paraminfo.Q;
   searchtype = paraminfo.localsearch;
   iteration = 1;
   
   ntours  = 0;
   localsearchtype = 2;
   [tourdistance, ~] = nntour(Distance, nE, ntours, localsearchtype,...
                              paraminfo.nnsize, paraminfo.nearneighbours );
   trail0 = 1.0 / (paraminfo.rho * tourdistance); 
   
   ETA = 1./(Distance + 0.1);
   TAU = ones(nE, nE) * trail0;   
   TOTAL = TAU.^alpha .* ETA.^beta;
   TOTAL(1:1+size(TOTAL,1):end) = trail0;
   
   bestant = [];
   besttour = inf;
   bestiteration = iteration;
   BESTANT0 = zeros(antpop, nE);
   BESTTOURS = zeros(antpop, 1);
   
   clusterinfo.nSize = paraminfo.nsize;       % number of nodes in cluster
   clusterinfo.nSectors = paraminfo.nsectors; % number of sectors
   clusterinfo.nE  = paraminfo.nE;
   clusterinfo.X = paraminfo.X;
   clusterinfo.Y = paraminfo.Y;
   clusterinfo.Distance = paraminfo.Distance;
   
   [CLUSTERS] = createclusters(clusterinfo);
   MAXCLUSTERS = max(CLUSTERS)';
   rhomax  = paraminfo.rho;
   rhomin = 0.01;
    
   while iteration <= maxiteration
       
       iteration = iteration + 1;
        
       % Construct the solution using the node clustering principle or solve the TSP
       [ANT0, TOURS] = constructtoursnc(CLUSTERS, TOTAL, MAXCLUSTERS, antpop, nE, Distance, nearneighbours, ETA, TAU, alpha, beta);
              
       if searchtype > 0
          [ANT0, TOURS] = localsearch(ANT0, TOURS, searchtype, antpop, Distance, nE, nnsize, nearneighbours);           
       end       
       ANT0 = ANT0';
       
       % update information
       [iterbesttour, minindex] = min(TOURS);
       iterbestant = ANT0(minindex,:);
        
       if iterbesttour < besttour
          besttour = iterbesttour;
          bestant = ANT0(minindex,:);
          bestiteration = iteration;         
        % Calculate the branching factor
       end  
       
       if (iteration - 1) == 1
            BESTANT0(:, 1:end) = ANT0(:, 1:end);
            BESTTOURS(:) = TOURS(:);
       else
           
           [maxtour, maxindex] = max(BESTTOURS);           
           if iterbesttour < maxtour && isempty(find(BESTTOURS == iterbesttour,1)) 
              BESTANT0(maxindex, :) = iterbestant(:);
              BESTTOURS(maxindex) = iterbesttour; 
           end
       end 
              
       % Global pheromone update using the best ant                
       ant = [bestant bestant(1)];
       % determine adaptive pheromone value 
       rho = adaptivepheromonevalue(ant, ANT0', nE, antpop, rhomax, rhomin);
              
       DTau = Q / besttour;
       for i = 1 : nE 
           % amount of pheromone added at path i to i+1 by the best ant
           j = ant(i); h = ant(i+1); 
         
           TAU(j, h) = (1 - rho) * TAU(j, h) + rho * DTau;
           TAU(h, j) = TAU(j, h);
           
           TOTAL(h, j) = (TAU(h, j)^alpha) * (ETA(h, j)^beta); 
           TOTAL(j, h) = TOTAL(h, j);
       end
       
       % Print some data at some point
       if maxiteration >= 200
           if mod(iteration, 250) == 0 
              fprintf('Best tour so far: %d. found at iteration %d.\n', besttour, bestiteration); 
           elseif mod(iteration, 50) ==0
              fprintf('ant system iteration: %d.\n', iteration);
           end
       else
           if mod(iteration, 10) == 0 
              fprintf('Best tour so far: %d. found at iteration %d.\n', besttour, bestiteration); 
           elseif mod(iteration, 5) ==0
              fprintf('ant system iteration: %d.\n', iteration);
           end
       end 
   end 
   
   res.bestant = bestant;
   res.besttour = besttour;
   res.bestiteration = bestiteration;
   res.TAU = TAU;
   res.TOTAL = TOTAL;
   res.ETA = ETA;  
   res.ants = BESTANT0;     
end

function newrho = adaptivepheromonevalue(ant, ANT0, nE, antpop, rhomax, rhomin)
   
   Hmin = -log2(1/nE);
   Hmax = -log2(1/(nE*antpop));  
   OCCURS = createlinkfrequencies(ANT0, nE, antpop);
   
   % Adaptive pheromone evaluation
   SUM = 0;   
   for i = 2 : nE           
       for j = 1 : i-1
           if  i ~= j
%              occur = occurrences(ant, j, i);
              occur = OCCURS(j,i);
              pij = occur/(nE * antpop);
              
              if pij ~= 0
                 SUM = SUM + pij * log2(pij);
              end
           end  
       end
       
       if -SUM < Hmax
           H = -SUM;
       else
           break;
       end 
   end
%   H = -SUM;
   newrho = abs( rhomin + (rhomax - rhomin) * ((H - Hmin)/(Hmax - Hmin)) );
end

function nooccur = occurrences(ant, i, j)
      
  nooccur = 0;
  to = ant(j);      
  while true           
     if to == ant(i)
        break; 
     end       
     i = i + 1;
     nooccur = nooccur + 1;  
  end
end

function [OCCURANCES] = createlinkfrequencies(ANT0, nE, antpop)

    OCCURANCES = zeros(nE,nE);
    for k = 1 : antpop      
       from = 1;
       to = 2;
       while to <= nE
           OCCURANCES(ANT0(from,k), ANT0(to,k)) = OCCURANCES(ANT0(from,k), ANT0(to,k)) + 1;
           OCCURANCES( ANT0(to,k), ANT0(from,k)) = OCCURANCES(ANT0(from,k), ANT0(to,k));
           from = to;
           to = to + 1;
       end
       OCCURANCES(ANT0(from,k), ANT0(1,k)) = OCCURANCES(ANT0(from,k), ANT0(1,k)) + 1;
       OCCURANCES(ANT0(1,k), ANT0(from,k)) = OCCURANCES(ANT0(from,k), ANT0(1,k)); 
   end 
end
