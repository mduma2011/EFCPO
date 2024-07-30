function [Dist] = funcfasttsp(ant_path,nE, VecDistance)
% funcfasttsp is a calculation for TSP
% VecDistance is a vector of the Distance matrix
% This allows for fast index in the case of larget TSPs
% Calculation of distance traveled by ant 'k'

    Dist = 0;
    for t = 1:nE
       row = ant_path(t); 
       col = ant_path(t + 1);
       Dist = Dist + VecDistance(col * nE - (nE - row)); 
    end
end

