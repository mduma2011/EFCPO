function R = randomnumbers(a,b, rows, cols)
% Find the random number between a nd b
% rows in the no of rows
% cols in the no of cols

R = a + (b-a).* rand(rows, cols);

end

