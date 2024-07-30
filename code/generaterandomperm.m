function randvec = generaterandomperm(N)
    
    randvec = (1:N)';
    totalassigned = 0;
    
    for i = 1:N
       rnd = rand;
       node = floor( rnd * (N - totalassigned));
       temp = randvec(i);
       randvec(i) = randvec(i + node);
       randvec(i + node) = temp;
       totalassigned = totalassigned + 1;
    end
end
