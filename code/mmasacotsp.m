function [res] = mmasacotsp(paraminfo)
   % Max-Min Ant System - Ant Colony Optimization
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
   
   iteration = 1;
   ntours = 0;
   
   % help create nearest neighbour tour distance                                  
   [tourdistance, ~] = nntour(Distance, nE, ntours, searchtype, nnsize, nearneighbours);
   TAUMAX  = 1.0 / (rho * tourdistance); 
   TAUMIN  = TAUMAX / (2 * nE); 
   TAU0 = TAUMAX;
   
   % Initialize pheromone trails
   ETA = 1./(Distance + 0.1);
   TAU = ones(nE, nE) * TAUMAX; 
  
   TOTAL = TAU.^alpha .* ETA.^beta;
   TOTAL(1:1+size(TOTAL,1):end) = TAU0;
   
   bestant = [];
   besttour = inf;
   bestiteration = iteration;
   BESTANT0 = zeros(antpop, nE);
   BESTTOURS = inf(antpop, 1);
   
   while iteration <= maxiteration
       
        iteration = iteration + 1;
           
         % Construct the solution or solve the TSP
        [ANT0, TOURS] = constructtours(TOTAL, nearneighbours, antpop, nnsize, nE, Distance);
              
         if searchtype > 0
          
            [ANT0, TOURS] = localsearch(ANT0, TOURS, searchtype, antpop, Distance, nE, nnsize, nearneighbours);
            
         end 
         ANT0 = ANT0';
       
         % update information
         [iterbesttour, minindex] = min(TOURS);
          iterbestant = ANT0(minindex,:);
         
         if iterbesttour < besttour
            besttour = iterbesttour;
            bestant = iterbestant;
            bestiteration = iteration; 
                                  
            if searchtype <= 0
               px = exp( log(0.05) / nE);
               den = floor((nnsize + 1) / 2);
               TAUMIN = 1 * (1 - px) / (px * den);  
               TAUMAX = 1 / (rho * besttour);
               TAUMIN = TAUMAX * TAUMIN;
            else
               TAUMAX = 1 / (rho * besttour); 
               TAUMIN = TAUMAX / (2 * nE);               
            end           
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
                  
         % Update pheromone information.evaporation
         if searchtype > 0
             
            for i = 1:nE      
              for j = 1:nnsize
                  temp = nearneighbours(i, j);  
                  TAU(i, temp) = (1 - rho) * TAU(i, temp);
                  TAU(i,j) = max(TAU(i,j), TAUMIN); 
              end
            end             
         else    
            TAU = (1 - rho) * TAU;
         end
         
         % Global pheromone update using iteration best ant                
         ant = [iterbestant iterbestant(1)];
         DTau = Q / iterbesttour;
         for i = 1 : nE 
            % Amount of pheromone added at path i to i+1 by best ant
            j = ant(i); h = ant(i+1);          
            TAU(j, h) = TAU(j, h) + DTau;
            TAU(h, j) = TAU(j, h);           
         end
                  
         % Check for pheromone trail limits for max -min algorithm
         for i = 1:nE      
            for j = 1:i
                TAU(i,j) = max(TAU(i,j), TAUMIN);
                TAU(i,j) = min(TAU(i,j), TAUMAX);                
                TAU(j,i) = TAU(i,j);               
            end      
         end    

         % Update the totals
         if searchtype > 0
             
            % Update the totals using nearest neigbhours nnsize              
            for i = 1:nE      
              for j = 1:nnsize
                  h = nearneighbours(i, j);                  
                  if TAU(i, h) < TAU(h, i)
                     TAU(h, i) = TAU(i, h); 
                  end    
                  
                  TOTAL(i, h) = TAU(i, h)^alpha * ETA(i, h)^beta;
                  TOTAL(h, i) = TOTAL(i, h) ;
              end
            end
           
         else
          
            % Update the totals using number of nodes nE  
            TOTAL = (TAU.^alpha) .* (ETA.^beta);
            TOTAL(1:1+size(TOTAL,1):end) = TAU0;
         end
                           
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
