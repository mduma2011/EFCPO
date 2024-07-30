function [vec1, vec2] = swap2(vec1, vec2, i, j )

   temp = vec1(i);
   vec1(i) = vec1(j);
   vec1(j) = temp;
   
   temp = vec2(i);
   vec2(i) = vec2(j);
   vec2(j) = temp;
end

