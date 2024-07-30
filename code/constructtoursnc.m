function [ANT0, TOURS] = constructtoursnc(CLUSTERS, TOTAL, MAXCLUSTERS, antpop, nE, Distance, nearneighbours, ETA, TAU, alpha, beta)
  %%----------------------------------------------------------------------
  %  Construct tours using a node clustering 
  %-----------------------------------------------------------------------
  [ANT0] = getantsolutionsnc(TOTAL, nearneighbours, CLUSTERS, MAXCLUSTERS, ETA, TAU, alpha, beta, antpop, nE);
            
  [TOURS] = computetours(ANT0, antpop, nE, Distance);

end

function [ANT0] = getantsolutionsnc(TOTAL, nearneighbours, CLUSTERS, MAXCLUSTERS, ETA, TAU, alpha, beta, antpop, nE)

    VISITED = false(nE, antpop); % mark the citites as not visited
    step = 1;                    % initally place all the ants on same node
    ANT0 = zeros(nE, antpop);
       
    for i = 1 : antpop           
       % node = randi([1 nE],1,1); 
        node = nearneighbours(randi([1 nE],1,1), 1);
        ANT0(step, i) = node;
        VISITED(node, i) = true;
    end 
        
    for step =  2 : nE            
        for k = 1 : antpop 
               
            % 1. visit the node
            lastnode = ANT0(step-1, k); 
            Clastnode = CLUSTERS(:,lastnode);
            nClusters = MAXCLUSTERS(lastnode);
            
            PROB = probabilitypercluster(k, Clastnode, nClusters, VISITED, ETA(:,lastnode), TAU(:,lastnode), alpha, beta, nE);
            
            [~, cluster] = max(PROB);
  
            % visit the next node from the list of nodes in the selected cluster
             candidates = find(Clastnode == cluster);
            
            if isempty(find(VISITED(candidates, k)== false,1)) && cluster <= nClusters
                error('Error: Something went wrong with the probality per cluster calculation');  
            end    
            
            candsize = length(candidates);            
            nextnode = visitnextnodenc(k, lastnode, VISITED, TOTAL,candidates, candsize);
            ANT0(step, k) = nextnode;
            VISITED(nextnode, k) = true; 
         end
    end    
end

function nextnode = visitnextnodenc(k, vk, VISITED, TOTAL,candidates, candsize)

   node = vk;
   PROB = zeros(candsize, 1); 
   SUM  = 0;
   for i = 1 : candsize

      candidatenode = candidates(i);                   
      if ~VISITED(candidatenode, k)  
          PROB(i) = TOTAL(node, candidatenode);
          SUM = SUM + PROB(i);
       end   
   end  
   
   if SUM <= 0
      error('Error: Sum is less than zero under visitnextnodenc'); 
   end    
          
   % This is a roulette selection if one fails to find a node
   crand  = rand;
   randomprob = crand;
   randomprob = randomprob * SUM;
   i = 1;
   partialSum = PROB(i);

   while partialSum <= randomprob && i < length(PROB)
      i = i + 1;
      partialSum = partialSum + PROB(i);  
   end
   temp = candidates(i);

   if VISITED(temp, k) 
       i = 1;
       while PROB(i) <= 0
          i = i + 1; 
       end    
       temp = candidates(i);
   end
   nextnode = temp;       
end

function PROB = probabilitypercluster(k, Clastnode, nClusters, VISITED, eta, tau, alpha, beta, nE)

     ETAc = zeros(nClusters, 1);
     TAUc = zeros(nClusters, 1);
     PROB = zeros(nClusters, 1);
     NC = ones(nClusters, 1);
     
     for node = 1 : nE
        
         if ~VISITED(node, k)
            
             cluster = Clastnode(node);
             ETAc(cluster) = ETAc(cluster) + eta(node);
             TAUc(cluster) = TAUc(cluster) + tau(node);
             NC(cluster) = NC(cluster) + 1;
         end    
     end    
     
     ETAc = ETAc ./ NC;
     TAUc = TAUc ./ NC;
     
     NUM = (ETAc.^alpha) .* (TAUc.^beta);
     DEN = sum(NUM);
      
     if DEN > 0
        PROB = NUM/DEN; 
     end
end

