function nextnode = visitnextnode(k, ANT0, VISITED, TOTAL, nearneighbours, step, nnsize, nE)
        
    node = ANT0(step-1, k);
    PROB = zeros(nnsize, 1); 
    SUM  = 0;
    for i = 1 : nnsize

      nneibhor = nearneighbours(node, i);                   
      if ~VISITED(nneibhor, k)  
          PROB(i) = TOTAL(node, nneibhor);
          SUM = SUM + PROB(i);
       end   
    end  
 
   % choose the next best
   if SUM <= 0
      
      nextcity = nE;
      node = ANT0(step-1, k);
      bestval = -inf;                   
      for city = 1 : nE

         if ~VISITED(city, k)
             if TOTAL(node, city) > bestval         
                nextcity = city;
                bestval =  TOTAL(node, city);
             end    
         end
      end     
      nextnode = nextcity;
      
   else  
       
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
       temp =  nearneighbours(node, i);
             
       if VISITED(temp, k) 
           i = 1;
           while PROB(i) <= 0
              i = i + 1; 
           end    
           temp =  nearneighbours(node, i);
       end
       nextnode = temp;       
   end 
end