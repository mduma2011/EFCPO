
function nextnode = visitnextneigbhournode(k, ANT0, VISITED, TOTAL, nearneighbours, step, nnsize, nE)
    
    currnode = ANT0(step-1, k);
    nextnode = nE + 1;
    bestvalue = -inf; 
    
    for i = 1 : nnsize
        
       tempnode = nearneighbours(currnode, i);
       if VISITED(tempnode, k) == false
           if TOTAL(currnode, tempnode) > bestvalue
              bestvalue = TOTAL(currnode, tempnode);
              nextnode = tempnode;
           end    
       end
    end    

   % choose the next best
   if nextnode == nE + 1

      nextcity = nE;
      currnode = ANT0(step-1, k);
      bestval = -inf;                   
      for city = 1 : nE

         if VISITED(city, k) == false                         
             if TOTAL(currnode, city) > bestval         
                nextcity = city;
                bestval =  TOTAL(currnode, city);
             end    
         end
      end      
      nextnode = nextcity;   
   end
end

