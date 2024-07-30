function [TOURS] = computetours(ANT0, antpop, nE, Distance)

     TOURS = zeros(antpop, 1);         
     for k = 1 : antpop
        
         ant = ANT0(:,k);
         TOURS(k) = functsp([ant; ant(1)],nE, Distance);
         
     end    
end