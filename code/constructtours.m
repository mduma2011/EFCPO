function [ANT0, TOURS] = constructtours(TOTAL, nearneighbours, antpop, nnsize, nE, Distance)
            
       [ANT0] = getantsolutions(TOTAL, nearneighbours, antpop, nnsize, nE);
  
       [TOURS] = computetours(ANT0, antpop, nE, Distance);
end

function [ANT0] = getantsolutions(TOTAL, nearneighbours, antpop, nnsize, nE)
          
       % One needs to avoid using nE in the column section
       % in case nE is too large which may affect the performance
       % of the matrix calculation or tour construction as a whole.
          
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
               
               % --- choose neighbour and move to the next node
               nextnode = visitnextnode(k, ANT0, VISITED, TOTAL, nearneighbours, step, nnsize, nE);
               ANT0(step, k) = nextnode;
               VISITED(nextnode, k) = true; 
           end
       end      
end