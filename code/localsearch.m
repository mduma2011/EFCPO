function [RANT0,RTOURLEN] = localsearch(ANT0, TOURLEN, type, antpop, Distance,nE, nnsize, nearneighbours)

  optinfo.Distance = Distance;
  optinfo.nE = nE;
  optinfo.nnsize = nnsize;
  optinfo.nearneighbours = nearneighbours;  
    
  for k = 1 : antpop

    optinfo.ant = ANT0(:,k);
    optinfo.tourlength = TOURLEN(k);
    switch type
    
        case 2
             
             [newant, noimprovements, newlen]  = twooptsearch(optinfo);             
             if noimprovements > 0
                ANT0(:,k) =  newant(:); 
                TOURLEN(k) = newlen;
             end   
        case 3
            
              [newant, noimprovements, newlen] = threeoptsearch(optinfo);              
              if noimprovements > 0
                 ANT0(:,k) =  newant(:); 
                 TOURLEN(k) = newlen;
              end   
        otherwise
           printf('local search not correctly specified, %d instead of 0, 2 or 3');           
    end    
  end
   
  RANT0 = ANT0;
  RTOURLEN = TOURLEN;
end

