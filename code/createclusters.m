function [CLUSTERS] = createclusters(clusterinfo)
%%-------------------------------------------------------
% Create a cluster of nodes by selecting nodes into
% clusters.
%--------------------------------------------------------

 nSize = clusterinfo.nSize;       % number of nodes in cluster
 nSectors = clusterinfo.nSectors; % number of sectors
 nE  = clusterinfo.nE;
 X = clusterinfo.X;
 Y = clusterinfo.Y;
 Distance = clusterinfo.Distance;
 CLUSTERS = zeros(nE,nE);
 VISITED  = false(nE,nE);
 
 if (1 <= nSize) && (nSize <= nE)
     
     for vi = 1 : nE
         
         VISITED(vi,vi) = true;
         cluster = 1;
                 
         for j = 1 : nSectors
             selectednode = findclosestnodeinsector(vi, VISITED, nSectors, j,nE,X,Y,Distance); 
             if selectednode > 0
                CLUSTERS(selectednode, vi) = cluster;
                VISITED(selectednode, vi) = true;
             end   
         end
         
         while ~isempty(find(VISITED(:, vi) == false,1))                     % theres no selected node
              while length(find(CLUSTERS(:, vi) == cluster))  < nSize
                                
                 selectednode = findclosestnode(vi, VISITED, nE, Distance);
                 if selectednode > 0
                    CLUSTERS(selectednode, vi) = cluster;
                    VISITED(selectednode, vi) = true;
                 else
                     break;
                 end         
              end
              cluster = cluster + 1;
         end         
     end 
 else
    error('nSize is less than 1 or greater that the number of nodes'); 
 end    
 
end

function closestnode = findclosestnode(vi, VISITED, nE, Distance)

%%------------------------------------------------------------------------
%  This one is easy, just find the closest node from the nodes connected to
%  vi, provided it's not already selected
%
%-------------------------------------------------------------------------

 mindist = inf;
 closestnode = -1;
 for node = 1 : nE
     
     if ~VISITED(node, vi)
        dist = Distance(node, vi);
        
        if mindist > dist
           mindist = dist;
           closestnode = node;
        end           
     end        
 end
 
end


function closestnode = findclosestnodeinsector(vi, VISITED, nSectors, j,nE,X,Y,Distance)
 
%%--------------------------------------------------------------------------------------
%  This implementation is from the paper "Adaptive Ant Colony Optimization 
%  with node clustering applied to the Travelling Salesman Problem by 
%  Petr Stodola et. al.". It's not clear how finding the closest is explained.
%  In fact, there is no explanation on how the sectors are created. It's
%  just a diagram, Figure 3, that shows you a figure with the sectors already 
%  created. 
%
%  I figured out that the best way to create the sectors is to split the 
%  sectors using angles and quadrants. In other words, let's the number
%  of sectos is 2, this can be used to split the sectors into 2 quagrants sepated by pi,  
%  i.e. 2 * pi / 2 =  pi. The 2 quadrants are separated by the angle pi.
%  The following formula,2 * pi /nSectors, can be used in the following way: 
%  1. for 2 sectors, the sectors are separated by 2*pi/2 =  pi
%  2. for 4 sectors, the sectors are separated by 2*pi/4 =  pi/2
%  3. for 6 sectors, the sectors are separated by 2*pi/6 =  pi/3
%  4. for 8 sectors, the sectors are separated by 2*pi/8 =  pi/4
%
%---------------------------------------------------------------------------------------

 rednode = vi;
 offsetx = X(rednode);
 offsety = Y(rednode);
 
 sectors = 0:2*pi/nSectors:2*pi;
 
 mindist = inf;
 closestnode = -1;
 
 % 1. iterate through each node. 
 for k = 1 : nE
     
    if ~VISITED(k, rednode)
        
       node = k;       
       % 2. We already know we working with sector j
       % check if the node is part of the sector
       x = X(node); 
       y = Y(node);
           
       dx = x - offsetx;
       dy = y - offsety;           
       if (x - offsetx) == 0 
           radval = 2*pi + asin(dy/abs(dy));  
       elseif (y - offsety) == 0
           radval = 2*pi + acos(dx/abs(dx)); 
       else
           radval = atan(dy/dx);              
           if radval < 0
              if dx < 0
                 radval = radval + pi;   
              elseif dy < 0    
                 radval = radval + 2*pi;
              end
           end   
       end   
           
      if sectors(j) < radval &&  radval < sectors(j+1)          
          % get the distance between rednode and node k
          dist = Distance(rednode, k);          
          % check if it is the the next closest node
          if mindist > dist
             mindist = dist;
             closestnode = k;
          end              
      end      
    end
 end
end

