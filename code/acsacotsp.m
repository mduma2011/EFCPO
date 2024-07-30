
function [res] = acsacotsp(paraminfo)
   % Ant System With Local Pheromone Update - Ant Colony Optimization
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
   q0 = paraminfo.q0;
   iteration = 1;
   
   ntours  = 0;
   [tourdistance, ~] = nntour(Distance, nE, ntours, searchtype,...
                              paraminfo.nnsize, paraminfo.nearneighbours );
   
   
   trail0  = 1.0 / (nE * tourdistance); 
      
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
           
       % Version 2. Construct the solution
       localrho = 0.1;       
       VISITED = false(nE, antpop); % mark the citites as not visited        
       step = 1;                    % initally place all the ants on same node
       ANT0 = zeros(nE, antpop);
       
       for i = 1 : antpop
           
          node = randi([1 nE],1,1);            
          ANT0(step, i) = node;
          VISITED(node, i) = true;           
       end 
       
       for step =  2 : nE            
           for k = 1 : antpop 
               
               % --- choose neighbour and move to the next node
               nextnode = selectnextnode(k, ANT0, VISITED, TOTAL, nearneighbours, step, nnsize, nE,q0);
               
               ANT0(step, k) = nextnode;
               VISITED(nextnode, k) = true;
               
               % --- local pheromone update              
               tour = ANT0(:, k);                   
               j = tour(step);
               h = tour(step - 1);
               TAU(h,j) = (1 - localrho) * TAU(h,j) + localrho * trail0;
               TAU(j,h) = TAU(h,j); 
               TOTAL(h,j) = TAU(h,j)^alpha * ETA(h, j)^beta;
               TOTAL(j,h) = TOTAL(h,j);                  
           end
       end
              
       % Version 2. Construct tours for each antk
        TOURS = zeros(antpop, 1);  
        step = nE;
        for k = 1 : antpop

             ant = ANT0(:, k);
             TOURS(k) = functsp([ant; ant(1)],nE, Distance);

             % --- local pheromone update
             tour = ANT0(:, k);                   
             j = tour(step);
             h = tour(step - 1);
             TAU(h, j) = (1 - localrho) * TAU(h, j) + localrho * trail0;
             TAU(j, h) = TAU(h, j); 
             TOTAL(h, j) = (TAU(h, j)^alpha) * (ETA(h, j)^beta);
             TOTAL(j, h) = TOTAL(h, j);           
       end   
                  
       if searchtype > 0          
          [ANT0,TOURS] = localsearch(ANT0, TOURS, searchtype, antpop, Distance, nE, nnsize, nearneighbours);           
       end    
       
       ANT0 = ANT0';
       
       % update information
       [iterbesttour, minindex] = min(TOURS);
         iterbestant = ANT0(minindex,:);
       
       if iterbesttour < besttour
          besttour = iterbesttour;
          bestant = ANT0(minindex,:);
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
            
       % Global pheromone update using the best ant                
       ant = [bestant bestant(1)];
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

function nextnode = selectnextnode(k, ANT0, VISITED, TOTAL, nearneighbours, step, nnsize, nE, q0)
     
   if rand < q0
       
       nextnode = visitnextneigbhournode(k, ANT0, VISITED, TOTAL, nearneighbours, step, nnsize, nE);       
       return;
   end
   
   nextnode = visitnextnode(k, ANT0, VISITED, TOTAL, nearneighbours, step, nnsize, nE);
end              


