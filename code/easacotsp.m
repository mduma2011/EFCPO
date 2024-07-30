
function [res] = easacotsp(paraminfo)
  % Elitist Ant System - Ant Colony Optimization
  % implementation

   alpha = paraminfo.alpha;
   beta = paraminfo.beta;
   rho = paraminfo.rho;
   nnsize = paraminfo.nnsize;
   nE = paraminfo.nE;
   Distance = paraminfo.Distance ;
   nearneighbours = paraminfo.nearneighbours;
   maxiteration = paraminfo.maxiteration;
   antpop = paraminfo.antpop;
   Q = paraminfo.Q;
   searchtype = paraminfo.localsearch;
   elitistants = paraminfo.elitistants;
   
   % nodesvisited = logical(zeros(antpop, nnsize));
   % phase = 0;
   iteration = 1;
   ntours  = 0;
   
   % help create nearest neighbour tour distance                                  
   [tourdistance, ~] = nntour(Distance, nE, ntours, searchtype, nnsize, nearneighbours);
   trail0  = 1.0 / (rho * tourdistance); 
   
   % Initialize pheromone trails
   ETA = 1./(Distance + 0.1);
   TAU = ones(nE, nE) * trail0; 
  
   TOTAL = TAU.^alpha .* ETA.^beta;
   TOTAL(1:1+size(TOTAL,1):end) = trail0;
   
   bestant = [];
   besttour = inf;
   bestiteration = iteration;
   BESTANT0 = zeros(antpop, nE);
   BESTTOURS = zeros(antpop, 1);
   
   while iteration <= maxiteration
       
         iteration = iteration + 1;
        
         % Construct the solution or solve the TSP
         [ANT0, TOURS] = constructtours(TOTAL, nearneighbours, antpop, nnsize, nE, Distance);

         if searchtype > 0
            [ANT0,TOURS] = localsearch(ANT0, TOURS, searchtype, antpop, Distance, nE, nnsize, nearneighbours);
         end
         ANT0 = ANT0';
       
         % update information
         [iterbesttour, minindex] = min(TOURS);
          iterbestant = ANT0(minindex,:);

         if iterbesttour < besttour
             besttour = iterbesttour;
             bestant = iterbestant;
             bestiteration = iteration; 
           %  Calculate the branching factor
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
           
          % update pheromone information. Evaporation
          TAU = (1 - rho) * TAU;
       
          % Global update pheromone per ant k
          for k = 1 : antpop
           
              ant = [ANT0(k,:) ANT0(k,1)];
              DTau = Q / TOURS(k);
              for i = 1 : nE 
                  % amount of pheromone added at path i to i+1 by ant k
                  j = ant(i); h = ant(i+1);
                  TAU(j, h) = TAU(j, h) + DTau;
                  TAU(h, j) = TAU(j, h);
              end
           end
           
           % Global update pheromone using the best ant
           weight = elitistants;
           ant = [bestant bestant(1)];
           DTau = weight / besttour;
           for i = 1 : nE                
               j = ant(i); h = ant(i+1);
               TAU(j, h) = TAU(j, h) + DTau;
               TAU(h, j) = TAU(j, h);                
           end  
            
           % Compute total information
           TOTAL = TAU.^alpha .* ETA.^beta;   
           TOTAL(1:1+size(TOTAL,1):end) = trail0;
           
           % Print some data at some point
           if maxiteration >= 200
              if mod(iteration, 500) == 0 
                 fprintf('Best tour so far: %d. found at iteration %d.\n', besttour, bestiteration); 
              elseif mod(iteration, 200) ==0
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