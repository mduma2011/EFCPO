function [newtour, noimprovements, newtourlen] = threeoptsearch(optinfo)

 nE  = optinfo.nE;
 Distance = optinfo.Distance;
 nnsize = optinfo.nnsize;
 nearneighbours = optinfo.nearneighbours;
 tour = [optinfo.ant; optinfo.ant(1)];
 newtourlen =  optinfo.tourlength;
    
 noimprovements = 0;
 N = nE;
 pos = zeros(N, 1);   
 dlb = false(N, 1);
 
 for i = 1 : N
     pos(tour(i)) = i;        
 end
    
 if ~isempty(find(pos == 0, 1))
      newtour = tour; 
      return;
 end    
    
 improvement = true;
 randvector = generaterandomperm(N); 
 between = false;  
 opt2flag =  false;
 cangotonextedge = false;
 
 while improvement == true
        
     movevalue = 0;  
     improvement = false;
      
     for L = 1 : N
         
         c1 = randvector(L);
         if dlb(c1) == true
            continue;  
         end
         
         opt2flag =  false;
         moveflag = 0;
         posc1 = pos(c1);
         succ1 = tour(posc1 + 1);
         
         % REDUNTANT CODE
%          if posc1 > 1 pc1 = tour(posc1 - 1);            
%          else  pc1 = tour(N); 
%          end    
         
         h = 0;                   % find the h-nearest neighbours
         while h < nnsize
            
             c2 = nearneighbours(c1, h + 1); % find the position of the second city
             posc2 = pos(c2);
             succ2 = tour(posc2 + 1);
             
             if posc2 > 1 
                prevc2 = tour(posc2 - 1);
             else
                prevc2 = tour(N); 
             end
                           
             diffs = 0;
             dissp = 0;
             radius  = Distance(c1, succ1);
             add1 = Distance(c1, c2);
             
             if radius > add1     % fixed radius neigbhour search is performed
                decreasebreaks = -radius - Distance(c2, succ2);
                diffs = decreasebreaks + add1 + Distance(succ1, succ2);
                diffp = -radius - Distance(c2, prevc2) +... 
                         Distance(c1, prevc2) + Distance(succ1, c2);
             else
                break; 
             end    
         
             if prevc2 == c1       % no exchange is possible  
                diffp = 0; 
             end    
             if (diffs < movevalue) || (diffp < movevalue) 
                
                 %improvement = true;
                 if diffs <= diffp
                    h1 = c1; h2 = succ1; h3 = c2; h4 = succ2;
                    movevalue = diffs;
                    opt2flag = true;
                    moveflag = 0;
                 else
                     h1 = c1; h2 = succ1; h3 = prevc2; h4 = c2;
                     movevalue = diffp;
                     opt2flag = true;
                     moveflag = 0;
                 end                
             end 
             
             g = 0;             % perform the innermost search
             while g < nnsize
                
                 c3 = nearneighbours(succ1, g + 1);
                 posc3 = pos(c3);
                 succ3 = tour(posc3 + 1);
                 
                 if posc3 > 1
                    prevc3 = tour(posc3 - 1); 
                 else
                    prevc3 = tour(N); 
                 end    
                 
                 if c3 == c1
                    g = g + 1;
                    continue;
                 else
                    add2 = Distance(succ1, c3);
                    
                    if (decreasebreaks + add1) < add2 % perform fixed radius search for innermost search
                       
                        if posc2 > posc1
                            
                           if posc3 <= posc2 && posc3 > posc1
                              between = true; 
                           else    
                              between = false;
                           end    
                        elseif posc2 < posc1  
                            
                            if posc3 > posc1 || posc3 < posc2
                               between = true; 
                            else    
                               between = false; 
                            end   
                        else    
                            fprintf("Error: pos1 %d == pos2 %d. This should not happen\n", posc1, posc2);
                        end
                        
                        
                        if between == true
                           % add edges (c1, c2), (c3, succ1), (pc3, succ2)
                           % for a valid tour.
                           gain = decreasebreaks - Distance(c3, prevc3) +... 
                                  add1 + add2 + Distance(prevc3, succ2);
                           
                           if gain < movevalue  %   check if there is an improvement
                              movevalue = gain;
                              opt2flag = false;
                              moveflag = 1;
                              
                              h1 = c1; h2 = succ1; h3 = c2; h4 = succ2; h5 = prevc3; h6 = c3;
                              [newtour, newpos, newdlb, newc1, newc2, newposc1, newposc2] = threeoptexchange(pos, tour, h1, h2, h3, h4, h5, h6,c1, c2, dlb,N, moveflag, opt2flag);
                              
                              tourlen = functsp(newtour, N, Distance);
                              if tourlen < newtourlen
                                 improvement = true; 
                                 newtourlen = tourlen; 
                                 tour(:) = newtour(:);
                                 pos(:) = newpos(:);
                                 dlb(:) = newdlb(:);
                                 c1 = newc1; c2 = newc2;
                                 posc1 = newposc1; posc2 = newposc2;
                                 noimprovements = noimprovements + 1;
                                 movevalue = 0;
                              end
                              
                              cangotonextedge = true;                             
                              break;
                           end 
                        else    
                            % Add to edges (c1,c2), (succ1, c3), (succ2,
                            % succ3)
                            gain = decreasebreaks - Distance(c3, succ3) +...
                                   add1 + add2 + Distance(succ2, succ3);
                               
                             if posc2 == posc3
                                gain = 20000;     % NOTE: This looks like a lazy workaround. Need to revisit this 
                             end    
                             
                             if gain < movevalue
                                movevalue = gain;
                                opt2flag = false;
                                moveflag = 2;
                               
                                h1 = c1; h2 = succ1; h3 = c2; h4 = succ2; h5 = c3; h6 = succ3;                               
                                [newtour, newpos, newdlb, newc1, newc2, newposc1, newposc2] = threeoptexchange(pos, tour, h1, h2, h3, h4, h5, h6,c1, c2, dlb,N, moveflag, opt2flag);
                                                              
                                tourlen = functsp(newtour, N, Distance);
                                if tourlen < newtourlen
                                   improvement = true; 
                                   newtourlen = tourlen; 
                                   tour(:) = newtour(:);
                                   pos(:) = newpos(:);
                                   dlb(:) = newdlb(:);
                                   c1 = newc1; c2 = newc2;
                                   posc1 = newposc1; posc2 = newposc2;
                                   noimprovements = noimprovements + 1;
                                   movevalue = 0;
                                end
                                
                                cangotonextedge = true;
                                break;
                             end    
                             
                             % or try to add edges (c1,c2), (succ1,c3), (pc2, pc3)
                             gain = -radius - Distance(prevc2, c2) - Distance(prevc3, c3) + ...
                                    add1 + add2 + Distance(prevc2, prevc3);
                                
                             if c3 == c2 || c2 == c1 || c1 == c3 || prevc2 == c1
                                 gain  = 2000000;    % Another lazy work around                                 
                             end    
                             
                             if gain < movevalue                                
                                movevalue = gain;
                                opt2flag = false;
                                moveflag = 3;
                                
                                h1 = c1; h2 = succ1; h3 = prevc2; h4 = c2; h5 = prevc3; h6 = c3;
                                [newtour, newpos, newdlb, newc1, newc2, newposc1, newposc2] = threeoptexchange(pos, tour, h1, h2, h3, h4, h5, h6,c1, c2, dlb,N, moveflag, opt2flag);

                                tourlen = functsp(newtour, N, Distance);
                                if tourlen < newtourlen
                                   improvement = true; 
                                   newtourlen = tourlen; 
                                   tour(:) = newtour(:);
                                   pos(:) = newpos(:);
                                   dlb(:) = newdlb(:);
                                   c1 = newc1; c2 = newc2;
                                   posc1 = newposc1; posc2 = newposc2;
                                   noimprovements = noimprovements + 1;
                                   movevalue = 0;
                                end
                                
                                cangotonextedge = true;
                                break;
                             end  
                             
                             % Remove edges (c1,succ1),(c2,pc2),(c3,succ3)
                             % and insert edges (c1,c2),(c3,succ1),(pc2,succ3)
                             
                             gain = -radius - Distance(prevc2, c2) - Distance(c3, succ3) + ...
                                    add1 + add2 + Distance(prevc2, succ3);
                                                     
                             if gain < movevalue      % check for improvements   
                                movevalue = gain;
                                opt2flag = false;
                                moveflag = 4;
                                
                                h1 = c1; h2 = succ1; h3 = prevc2; h4 = c2; h5 = c3; h6 = succ3;
                                [newtour, newpos, newdlb, newc1, newc2, newposc1, newposc2] = threeoptexchange(pos, tour, h1, h2, h3, h4, h5, h6,c1, c2, dlb,N, moveflag, opt2flag);
 
                                tourlen = functsp(newtour, N, Distance);
                                if tourlen < newtourlen
                                   improvement = true; 
                                   newtourlen = tourlen; 
                                   tour(:) = newtour(:);
                                   pos(:) = newpos(:);
                                   dlb(:) = newdlb(:);
                                   c1 = newc1; c2 = newc2;
                                   posc1 = newposc1; posc2 = newposc2;
                                   noimprovements = noimprovements + 1;
                                   movevalue = 0;
                                end
                                
                                cangotonextedge = true;
                                break;
                             end
                        end
                        
                    else
                         g = nnsize + 1;
                    end 
                     
                 end 
                 g = g + 1;
             end
             
             if cangotonextedge
                cangotonextedge = false; 
                break;   
             end  
             h = h + 1;             
         end
         
         if ~(moveflag > 0 || opt2flag == true)           
            dlb(c1) = true; 
         end    
     end 
 end 
 
 newtour = tour(1:end-1);
 newtourlen = functsp([newtour; newtour(1)], N, Distance);
end

function [newtour, newpos, newdlb, c1, c2, posc1, posc2] = threeoptexchange(pos, tour, h1, h2, h3, h4, h5, h6, oldc1, oldc2, dlb, N, moveflag, opt2flag)

 val = zeros(1, 3); 
 htour = ones(1, N);
 hhtour = ones(1, N); 
 c1 = oldc1;
 c2 = oldc2;
 posc1 = pos(c1);
 posc2 = pos(c2);
 
% Perfom the 3-opt exchange, otherwise do the 2opt exchange
if moveflag > 0
    
   dlb(h1) = false; dlb(h2) = false; dlb(h3) = false;
   dlb(h4) = false; dlb(h5) = false; dlb(h6) = false;
   posc1 = pos(h1); posc2 = pos(h3); posc3 = pos(h5);

   if moveflag == 4

       if posc2 > posc1
          n1 = posc2 - posc1; 
       else
          n1 = N - (posc1 - posc2); 
       end   
       if posc3 > posc2
          n2 = posc3 - posc2; 
       else
          n2 = N - (posc2 - posc3); 
       end
       if posc1 > posc3
          n3 = posc1 - posc3; 
       else
          n3 = N - (posc3 - posc1); 
       end

       % n1: len h2-h3, n2: len h4-h5, n3: len h6-h1
       val(1) = n1; val(2) = n2; val(3) = n3;

       % Order the partial tours
       h = 0;
       help = -Inf;
       for gi = 1 : 3
           if help < val(gi)
              help = val(gi);
              h = gi;
           end                   
       end

       % order partial tours based on length
       if h == 0

          % copy part from pos(h4) to pos(h5)
          j = pos(h4);
          h = pos(h5);
          i = 1;
          htour(i) = tour(j);
          n1 = 1;
          while j ~= h
             i = i + 1;
             j = j + 1;
             if j >= N + 1
                j = 1;  
             end    
             htour(i) = tour(j); 
             n1 = n1 + 1;
          end    

          % 1- copy partial tour 3 in new position
          j = pos(h4);
          i = pos(h6);
          tour(j) = tour(i); 
          pos(tour(i)) = j;
          while i ~= posc1
              i = i + 1;
              if i >= N + 1
                 i = 1; 
              end 
              j = j + 1;
              if j >= N + 1
                 j = 1; 
              end 
              tour(j) = tour(i);
              pos(tour(i)) = j;
          end 

          % 2 - Copy stored part from htour
          j = j + 1;
          if j >= N + 1
             j = 1; 
          end 
          for i = 1:n1
             tour(j) = htour(i);
             pos(htour(i)) = j;
             j = j + 1;
             if j >= N + 1
                j = 1; 
             end                          
          end
          tour(N+1) = tour(1);

       elseif h == 1   

          % copy part from pos(h6) to pos(h1)
          j = pos(h6);
          h = pos(h1);
          i = 1;
          htour(i) = tour(j);
          n1 = 1;
          while j ~= h
             i = i + 1;
             j = j + 1;
             if j >= N + 1
                j = 1;  
             end    
             htour(i) = tour(j); 
             n1 = n1 + 1;
          end    

          % 1- copy partial tour 3 in new position
          j = pos(h6);
          i = pos(h2);
          tour(j) = tour(i); 
          pos(tour(i)) = j;
          while i ~= posc2
              i = i + 1;
              if i >= N + 1
                 i = 1; 
              end 
              j = j + 1;
              if j >= N + 1
                 j = 1; 
              end 
              tour(j) = tour(i);
              pos(tour(i)) = j;
          end 

          % 2 - Copy stored part from htour
          j = j + 1;
          if j >= N + 1
             j = 1; 
          end 
          for i = 1:n1
             tour(j) = htour(i);
             pos(htour(i)) = j;
             j = j + 1;
             if j >= N + 1
                j = 1; 
             end                          
          end
          tour(N+1) = tour(1);

       elseif h == 2 

          % copy part from pos(h2) to pos(h3)
          j = pos(h2);
          h = pos(h3);
          i = 1; 
          htour(i) = tour(j);
          n1 = 1;
          while j ~= h
             i = i + 1;
             j = j + 1;
             if j >= N + 1
                j = 1;  
             end    
             htour(i) = tour(j); 
             n1 = n1 + 1;
          end    

          % 1- copy partial tour 3 in new position
          j = pos(h2);
          i = pos(h4);
          tour(j) = tour(i); 
          pos(tour(i)) = j;
          while i ~= posc3
              i = i + 1;
              if i >= N + 1
                 i = 1; 
              end 
              j = j + 1;
              if j >= N + 1
                 j = 1; 
              end 
              tour(j) = tour(i);
              pos(tour(i)) = j;
          end 

          % 2 - Copy stored part from htour
          j = j + 1;
          if j >= N + 1
             j = 1; 
          end 
          for i = 1:n1
             tour(j) = htour(i);
             pos(htour(i)) = j;
             j = j + 1;
             if j >= N + 1
                j = 1; 
             end                          
          end
          tour(N+1) = tour(1);
       end

   elseif moveflag == 1

       if posc2 > posc3  
          n1 = posc2 - posc3; 
       else
          n1 = N - (posc3 - posc2); 
       end   
       if posc3 > posc1
          n2 = posc3 - posc1 + 1; 
       else
          n2 = N - (posc1 - posc3 + 1); 
       end
       if posc2 > posc1
          n3 = N - (posc2 - posc1 + 1);
       else
          n3 = posc1 - posc2 + 1; 
       end

       % n1: len h6-h3, n2: len h5-h2, n3: len h1-h3
       val(1) = n1; val(2) = n2; val(3) = n3;

       % Order the partial tours
       h = 0;
       help = -Inf;
       for gi = 1 : 3                   
           if help < val(gi)
              help = val(gi);
              h = gi;
           end                   
       end

       % order partial tours based on length
       if h == 0

          % copy part from pos(h4) to pos(h5)
          j = pos(h5);
          h = pos(h2);
          i = 1;
          htour(i) = tour(j);
          n1 = 1;
          while j ~= h
             i = i + 1;
             j = j - 1;
             if j < 1
                j = N;  
             end    
             htour(i) = tour(j); 
             n1 = n1 + 1;
          end    

          % 1- copy partial tour 3 in new position
          j = pos(h1);
          h = pos(h4);
          i = 1;
          hhtour(i) = tour(j); 
          n2 = 1;
          while j ~= h
              i = i + 1;
              j = j - 1;
              if j < 1
                 j = N; 
              end 
              hhtour(i) = tour(j);
              n2 = n2 + 1;
          end 

          j = pos(h4);
          for i = 1:n2
             tour(j) = hhtour(i);
             pos(hhtour(i)) = j;
             j = j + 1;
             if j >= N +1
                j = 1; 
             end    
          end    

          for i = 1:n1
             tour(j) = htour(i);
             pos(htour(i)) = j;
             j = j + 1;
             if j >= N + 1
                j = 1; 
             end                          
          end
          tour(N+1) = tour(1);

       elseif h == 1   

          % copy part from pos(h3) to pos(h6)
          j = pos(h3);
          h = pos(h6);
          i = 1;
          htour(i) = tour(j);
          n1 = 1;
          while j ~= h
             i = i + 1;
             j = j - 1;
             if j < 1
                j = N;  
             end    
             htour(i) = tour(j); 
             n1 = n1 + 1;
          end    

          j = pos(h6);
          i = pos(h4);

          tour(j) = tour(i); 
          pos(tour(i)) = j;
          while i ~= posc1
              i = i + 1;
              j = j + 1;
             if j >= N +1
                j = 1; 
             end 
             if i >= N +1
                i = 1; 
             end
             tour(j) = tour(i);
             pos(tour(i)) = j;
          end 

          % Copy stored part from the htour
          j = j + 1;
          if j >= N + 1
             j = 1; 
          end    
          i = 1;
          tour(j) = htour(i);
          pos(htour(i)) = j; 
          while j ~= posc1                          
             j = j + 1;
             if j >= N +1
                j = 1; 
             end 
             i = i + 1;
             tour(j) = htour(i);
             pos(htour(i)) = j;
          end                      
          tour(N+1) = tour(1);

       elseif h == 2 
          %Check done.
          % copy part from pos(h2) to pos(h5)
          j = pos(h2);
          h = pos(h5);
          i = 1;
          htour(i) = tour(j);
          n1 = 1;
          while j ~= h
             i = i + 1;
             j = j + 1;
             if j >= N + 1
                j = 1;  
             end    
             htour(i) = tour(j); 
             n1 = n1 + 1;
          end    

          j = posc2;
          h = pos(h6);
          i = 1;
          hhtour(i) = tour(j);
          n2 = 1;
          while j ~= h
             i = i + 1;
             j = j - 1;
             if j < 1
                j = N;  
             end    
             hhtour(i) = tour(j); 
             n2 = n2 + 1;
          end

          j = pos(h2);
          for i = 1:n2
              tour(j) = hhtour(i);
              pos(hhtour(i)) = j;
              j = j + 1;
              if j >= N +1
                 j = 1; 
              end    
          end    

          % copy the stored part from htour
          for i = 1:n1
              tour(j) = htour(i);
              pos(htour(i)) = j;
              j = j + 1;
              if j >= N +1
                 j = 1; 
              end    
          end    
          tour(N + 1) = tour(1);
       end 

   elseif moveflag == 2

       if posc1 > posc3  
          n1 = posc1 - posc3; 
       else
          n1 = N - (posc3 - posc1); 
       end   
       if posc3 > posc2
          n2 = posc3 - posc2; 
       else
          n2 = N - (posc2 - posc3); 
       end
       if posc2 > posc1
          n3 = posc2 - posc1;
       else                       
          n3 = N - (posc1 - posc2); 
       end

       % n1: len h6-h3, n2: len h5-h2, n3: len h1-h3
       val(1) = n1; val(2) = n2; val(3) = n3;

       % Order the partial tours
       h = 0;
       help = -Inf;
       for gi = 1 : 3                   
           if help < val(gi)
              help = val(gi);
              h = gi;
           end                   
       end

       % order partial tours based on length
       if h == 0

          % copy part from pos(h4) to pos(h5)
          j = pos(h3);
          h = pos(h2);
          i = 1;
          htour(i) = tour(j);
          n1 = 1;
          while j ~= h
             i = i + 1;
             j = j - 1;
             if j < 1
                j = N;  
             end    
             htour(i) = tour(j); 
             n1 = n1 + 1;
          end    

          j = pos(h5);
          h = pos(h4);
          i = 1;
          hhtour(i) = tour(j); 
          n2 = 1;
          while j ~= h
              i = i + 1;
              j = j - 1;
              if j < 1
                 j = N; 
              end 
              hhtour(i) = tour(j);
              n2 = n2 + 1;
          end 

          j = pos(h2);
          for i = 1:n1
             tour(j) = htour(i);
             pos(htour(i)) = j;
             j = j + 1;
             if j >= N+1
                j = 1; 
             end    
          end    

          for i = 1:n2
             tour(j) = hhtour(i);
             pos(hhtour(i)) = j;
             j = j + 1;
             if j >= N + 1
                j = 1; 
             end                          
          end
          tour(N+1) = tour(1);

       elseif h == 1   

          % copy part from pos(h2) to pos(h3)
          j = pos(h2);
          h = pos(h3);
          i = 1;
          htour(i) = tour(j);
          n1 = 1;
          while j ~= h
             i = i + 1;
             j = j + 1;
             if j >= N + 1
                j = 1;  
             end    
             htour(i) = tour(j); 
             n1 = n1 + 1;
          end    

          j = pos(h1);
          h = pos(h6);
          i = 1;
          hhtour(i) = tour(j); 
          n2 = 1;                      
          while j ~= h
              i = i + 1;
              j = j - 1;
              if j < 1
                 j = N; 
              end 
              hhtour(i) = tour(j);
              n2 = n2 + 1;
          end 

          j = pos(h6);
          for i = 1:n1
             tour(j) = htour(i);
             pos(htour(i)) = j;
             j = j + 1;
             if j >= N+1
                j = 1; 
             end    
          end    

          for i = 1:n2
             tour(j) = hhtour(i);
             pos(hhtour(i)) = j;
             j = j + 1;
             if j >= N + 1
                j = 1; 
             end                          
          end           
          tour(N+1) = tour(1);

       elseif h == 2 
          %Check done
          % copy part from pos(h1) to pos(h6)
          j = pos(h1);
          h = pos(h6);
          i = 1;
          htour(i) = tour(j);
          n1 = 1;
          while j ~= h
             i = i + 1;
             j = j - 1;
             if j < 1
                j = N;  
             end    
             htour(i) = tour(j); 
             n1 = n1 + 1;
          end    

          j = pos(h4);
          h = pos(h5);
          i = 1;
          hhtour(i) = tour(j);
          n2 = 1;
          while j ~= h
             i = i + 1;
             j = j + 1;
             if j >= N + 1
                j = 1;  
             end    
             hhtour(i) = tour(j); 
             n2 = n2 + 1;
          end

          j = pos(h4);                      
          for i = 1:n1
              tour(j) = htour(i);
              pos(htour(i)) = j;
              j = j + 1;
              if j >= N +1
                 j = 1; 
              end    
          end

          % copy stored part from hhtour
          for i = 1:n2
              tour(j) = hhtour(i);
              pos(hhtour(i)) = j;
              j = j + 1;
              if j >= N +1
                 j = 1; 
              end    
          end                       
          tour(N + 1) = tour(1);
       end               

   elseif moveflag == 3 

       if posc1 > posc3  
          n1 = posc1 - posc3; 
       else
          n1 = N - (posc3 - posc1); 
       end   
       if posc3 > posc2
          n2 = posc3 - posc2; 
       else
          n2 = N - (posc2 - posc3); 
       end
       if posc2 > posc1
          n3 = posc2 - posc1;
       else                       
          n3 = N - (posc1 - posc2); 
       end

       val(1) = n1; val(2) = n2; val(3) = n3;

       % Order the partial tours
       h = 0;
       help = -Inf;
       for gi = 1 : 3                   
           if help < val(gi)
              help = val(gi);
              h = gi;
           end                   
       end

       % order partial tours based on length
       if h == 0

          % copy part from pos(h2) to pos(h3)
          j = pos(h3);
          h = pos(h2);
          i = 1;
          htour(i) = tour(j);
          n1 = 1;
          while j ~= h
             i = i + 1;
             j = j - 1;
             if j < 1
                j = N;  
             end    
             htour(i) = tour(j); 
             n1 = n1 + 1;
          end    

          j = pos(h2);
          h = pos(h5);
          i = pos(h4);
          tour(j) = h4; 
          pos(h4) = j;
          while i ~= h
              i = i + 1;
              if i >= N + 1
                 i = 1; 
              end
              j = j + 1;
              if j >= N + 1
                 j = 1; 
              end 
              tour(j) = tour(i);
              pos(tour(i)) = j;
          end 

          j = j + 1;
          if j >= N + 1
             j = 1;
          end                      
          for i = 1:n1
             tour(j) = htour(i);
             pos(htour(i)) = j;
             j = j + 1;
             if j >= N+1
                j = 1; 
             end    
          end    
          tour(N+1) = tour(1);

       elseif h == 1   

          % copy part from pos(h3) to pos(h2)
          j = pos(h3);
          h = pos(h2);
          i = 1;
          htour(i) = tour(j);
          n1 = 1;
          while j ~= h
             i = i + 1;
             j = j - 1;
             if j < 1
                j = N;  
             end    
             htour(i) = tour(j); 
             n1 = n1 + 1;
          end    

          j = pos(h6);
          h = pos(h1);
          i = 1;
          hhtour(i) = tour(j); 
          n2 = 1;                      
          while j ~= h
              i = i + 1;
              j = j + 1;
              if j >= N + 1
                 j = 1; 
              end 
              hhtour(i) = tour(j);
              n2 = n2 + 1;
          end 

          j = pos(h6);
          for i = 1:n1
             tour(j) = htour(i);
             pos(htour(i)) = j;
             j = j + 1;
             if j >= N+1
                j = 1; 
             end    
          end    

          for i = 1:n2
             tour(j) = hhtour(i);
             pos(hhtour(i)) = j;
             j = j + 1;
             if j >= N + 1
                j = 1; 
             end                          
          end           
          tour(N+1) = tour(1);

       elseif h == 2 

          % copy part from pos(h4) to pos(h5)
          j = pos(h5);
          h = pos(h4);
          i = 1;
          htour(i) = tour(j);
          n1 = 1;
          while j ~= h
             i = i + 1;
             j = j - 1;
             if j < 1
                j = N;  
             end    
             htour(i) = tour(j); 
             n1 = n1 + 1;
          end    

          j = pos(h1);
          h = pos(h6);
          i = 1;
          hhtour(i) = tour(j);
          n2 = 1;
          while j ~= h
             i = i + 1;
             j = j - 1;
             if j < 1
                j = N;  
             end    
             hhtour(i) = tour(j); 
             n2 = n2 + 1;
          end

          j = pos(h4); 
          % copy stored part from htour
          for i = 1:n1
              tour(j) = htour(i);
              pos(htour(i)) = j;
              j = j + 1;
              if j >= N +1
                 j = 1; 
              end    
          end

          % copy stored part from hhtour
          for i = 1:n2
              tour(j) = hhtour(i);
              pos(hhtour(i)) = j;
              j = j + 1;
              if j >= N +1
                 j = 1; 
              end    
          end                       
          tour(N + 1) = tour(1);
       end      
   else
       fprintf("Error: move flag not set correctly");
       return; % exit
   end
end  

if opt2flag == true

    dlb(h1) = false; dlb(h2) = false; 
    dlb(h3) = false; dlb(h4) = false;

    if pos(h3) < pos(h1)
       temp = h1; h1 = h3; h3 = temp;                     
       temp = h2; h2 = h4; h4 = temp; 
    end

    M = floor((N / 2) + 1);
    if (pos(h3)-pos(h2)) < M
        % reverse inner portion from pos(h2) to pos(h3)
        i = pos(h2); j = pos(h3);
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
        % reverse outer portion from pos(h4) to pos(h1)
        i = pos(h1); j = pos(h4);
        if j > i
           temp = N - (j - i) + 1; 
        else
           temp = (i - j) + 1; 
        end                    
        temp = temp / 2;

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
               i = N; 
            end    
            if j >= N + 1
               j = 1; 
            end    
        end
        tour(N+1) = tour(1);
    end                
 end
   
 newtour = tour;
 newpos  = pos;
 newdlb = dlb;
 posc1 = pos(c1); 
 posc2 = pos(c2);
end