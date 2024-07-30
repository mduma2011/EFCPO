function [vec1, vec2] = quicksort(vec1, vec2, left, right )
 
  if left >= right
     return; 
  end    
 
  [vec1, vec2] = swap2(vec1, vec2, left, floor((left + right)/2));
  last = left;
  
  for k = left+1:right
     if vec1(k) < vec1(left)
        last = last + 1; 
        [vec1, vec2] = swap2(vec1, vec2, last, k);  
     end
  end
  [vec1, vec2] = swap2(vec1, vec2, left, last);
  [vec1, vec2] = quicksort(vec1, vec2, left, last );
  [vec1, vec2] = quicksort(vec1, vec2, last + 1, right );
 
end

