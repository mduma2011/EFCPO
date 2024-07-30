function [tourdistance, ntours] = nntour(Distance, nE, ntours, searchtype, nnsize, nearneighbours)
   % This function helps create a tour length for a
   % random ant. This is to help initialise some of the 
   % parameters of the aco algorithm. 
   
   nodesvisited = false(nE,1);
   ant = zeros(nE, 1);
     
   nodei = randi([1 nE],1,1);  % place an ant at the start 
   ant(1) = nodei;
   tours = inf;
   nodesvisited(nodei) = true;
   
   for phase = 2:nE
      
       nextnode = nE;
       currnode = ant(phase-1);
       mindistance = inf;
       
       for node = 1:nE 
           
           if ~nodesvisited(node) 
               
               if Distance(currnode, node) < mindistance
                  mindistance =  Distance(currnode, node);
                  nextnode = node;
               end    
           end    
       end
       ant(phase) = nextnode;       
       nodesvisited(nextnode) = true;       
   end
      
   if searchtype > 0
         
       [~,tours] = localsearch(ant, tours, searchtype, 1, Distance, nE, nnsize, nearneighbours);
       tourdistance =  tours;
   else
       tourdistance = functsp([ant ant(1)],nE, Distance);     
   end   
   ntours  = ntours + 1;   
end
