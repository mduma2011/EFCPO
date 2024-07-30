function [loc] = roulette(probfitness)
   %ROULETTE wheel selection criteria
   partsum = 0;
   j = 1;
   sumfit = sum(probfitness);
   rn = rand * sumfit;
   
   while j <= length(probfitness)
      partsum = partsum + probfitness(j);
      if partsum >= rn
         break; 
      end
      j = j + 1;
   end
   
   if j > length(probfitness)
     j = length(probfitness);  
   end    
   loc = j;   
end