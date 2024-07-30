function [nearestnodeslist] = nnlist(Distance, nE, nnsize)
%  Distance: Nearest neighboour list
%  nE: number of ants 
  
  nearestnodeslist = zeros(nE, nnsize); 
  for node = 1:nE  
    
     % copy distances for nodes to a distance-vec 
     distancevec = Distance(node, :);
     nodevec = 1:nE;        
 
     distancevec(node) = Inf;
     [~, nodevec] = quicksort(distancevec, nodevec,1, nE);
     
     nearestnodeslist(node,:) = nodevec(1:nnsize);
  end
 
end

