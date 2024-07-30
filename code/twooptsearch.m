function [newtour, noimprovements, newtourlen] = twooptsearch(optinfo)

    nE = optinfo.nE;
    Distance = optinfo.Distance;
    nnsize = optinfo.nnsize;
    nearneighbours = optinfo.nearneighbours;
    tour = [optinfo.ant; optinfo.ant(1)];
       
    noimprovements = 0;
    noexchanges = 0;    
    N = nE;
    pos = zeros(N, 1);   
    dlb = false(N, 1);
    
    for i = 1 : N
        pos(tour(i)) = i;        
    end
        
    improvement = true;
    randvector = generaterandomperm(N); 
                      
    while improvement
        
        improvement = false;
        for L = 1 : N
            
            exchanged = false;
            c1 = randvector(L);            
            if dlb(c1)
               continue;  
            end
            
            posc1 = pos(c1);
            succ1 = tour(posc1 + 1);
            radius = Distance(c1, succ1);
            
            % search for posc1's nearest neigbhour
            
            for h = 1 : nnsize
                c2 = nearneighbours(c1, h);
                if radius > Distance(c1, c2)
                    
                   succ2 = tour(pos(c2) + 1);                  
                   gain = -radius + Distance(c1, c2) +...
                          Distance(succ1, succ2) - Distance(c2, succ2);
                   if gain < 0
                      h1 = c1; h2 = succ1;
                      h3 = c2; h4 = succ2;
                      
                      improvement = true;
                      noexchanges = noexchanges + 1;
                      dlb(h1) = false; dlb(h2) = false;
                      dlb(h3) = false; dlb(h4) = false;       
                      [newtour, newpos, ~, ~] = twooptexchange(pos, tour, h1, h2, h3, h4); 
                      
                      tour(:) = newtour(:);
                      pos(:) = newpos(:);  
                      exchanged = true;
                      break;                   
                   end                    
                else
                    break;
                end
            end
            
            if exchanged 
               continue; 
            end    
            
            % find the next c1's h-nearestt neighbours, use prevc1
            if posc1 > 1
               prevc1 = tour(posc1 - 1); 
            else
               prevc1 = tour(N); 
            end
            
            radius = Distance(prevc1, c1);
                        
            for h = 1 : nnsize
               c2 = nearneighbours(c1, h);
               if radius > Distance(c1, c2)
                  posc2 = pos(c2);
                  if posc2 > 1
                     prevc2 = tour(posc2 - 1); 
                  else
                     prevc2 = tour(N); 
                  end 
                  
                  if prevc2 == c1
                     continue;
                  end   
                  if prevc1 == c2
                     continue;
                  end
                                   
                  gain = -radius + Distance(c1, c2)+...
                          Distance(prevc1, prevc2) - Distance(prevc2, c2);
                  if gain < 0
                     h1 = prevc1; h2 = c1;
                     h3 = prevc2; h4 = c2;
                                        
                     improvement = true;
                     noexchanges = noexchanges + 1;
                     dlb(h1) = false; dlb(h2) = false;
                     dlb(h3) = false; dlb(h4) = false;       
                     [newtour, newpos, ~, ~] = twooptexchange(pos, tour, h1, h2, h3, h4); 

                     tour(:) = newtour(:);
                     pos(:) = newpos(:);
                     exchanged = true;
                     break;                    
                  end                  
               else
                  break; 
               end                
            end    
           
           if exchanged 
               continue; 
           end    
            
           dlb(c1) = true;
        end 
        
        if improvement
           noimprovements =  noimprovements + 1; 
        end       
    end
    
    newtour = tour(1:end-1);
    newtourlen = functsp([newtour; newtour(1)], N, Distance);
end

function [newtour, newpos, c1, c2] = twooptexchange(pos, tour, h1, h2, h3, h4)

     len = length(tour); 
     N = len;
     if pos(h3) < pos(h1)
        % perform swap 
        temp = h1; h1 = h3; h3 = temp;
        temp = h2; h2 = h4; h4 = temp;
     end   
        
     M = floor((N /2) + 1);
     
     if (pos(h3) - pos(h2)) < M
         % reverse inner part
         % 1. take route[0] to route[i] and add them in normal order to newtour.
         % 2. take route[i+1] to route[j] and add them in reverse order to newtour.
         % 3. take route[j+1] to route[end] and add them to in normal order to newtour.
         i = pos(h2); 
         j = pos(h3);

         while i < j
             c1 = tour(i);
             c2 = tour(j);
             tour(i) = c2;
             tour(j) = c1;
             pos(c1) = j;
             pos(c2) = i;
             i = i + 1;
             j = j - 1;     
         end
     else
         % reverse outer part
         % 1. take route[0] to route[i] and add them in reverse order to newtour.
         % 2. take route[i+1] to route[j] and add them in normal order to newtour.
         % 3. take route[j+1] to route[end] and add them to in reverse order to newtour.
         i = pos(h1); 
         j = pos(h4);
         if j > i 
            temp = N - (j - i) + 1; 
         else
            temp = (i - j) + 1; 
         end    
         
         temp = floor(temp / 2);
         
         for h = 1 : temp   
             c1 = tour(i);
             c2 = tour(j);
             tour(i) = c2;
             tour(j) = c1;
             pos(c1) = j;
             pos(c2) = i;  
                                                     
             i = i - 1;
             j = j + 1;
                  
             if i < 1
                i = len-1; 
             end             
             if j >= len 
                j = 1;            
             end
         end  
         tour(N) = tour(1);                
     end
       
     newtour = tour; 
     newpos = pos;
end
