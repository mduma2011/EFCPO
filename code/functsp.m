function [Dist] = functsp(ant_path,nE, Distance)
% functsp Distance calculation for TSP
% Calculation of distance traveled by ant 'k'

    Dist = 0;
    for t = 1:nE
       Dist = Dist + Distance(ant_path(t), ant_path(t + 1)); 
    end
end

