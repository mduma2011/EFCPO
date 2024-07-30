function [newalphas, newbetas, logevidences, errors, respercycle] = evidenceacotsp(evinfo, paraminfo)

% The Evidence Framework Driven Control Parameter Optimisation Algorithm (EFCPO)
% This is an implentaion done by M. Duma 2024.
% 
% In this implementation, rather than training a neural network (NN) and
% checking if that that NN model is good or not, or deriving
% hyperparameters alpha and beta for that NN model, we derive a hyperparameter 
% for an ACO model.
% 
% To check if it works maybe look at 2 ants buiding a solution per ACO and
% compare their log evidence value.

% We train the whole ACO model and the adujst its hypeparameter as we go. 
% this method calls the ACO model as an outer algorithm rather than the ACO 
% invoking it.

% We also retrain the model using the previous pheromone values TAU and
% TOTAL.
%
%--------------------------------------------------------------------------
% Initialisation of variables
%--------------------------------------------------------------------------

if(isempty(evinfo.alphas))
   error('Alphas should be populated.'); 
end

if(evinfo.cycles <= 0)
   error('Cycles cannot be zero or negative.'); 
end

cycles = evinfo.cycles;
ants = evinfo.ants;
X = evinfo.X;
Y = evinfo.Y;
D = evinfo.Distance; 
ETA = evinfo.ETA;
nE = evinfo.nE;
Error = evinfo.Error;
benchmark = evinfo.benchmark;
weightlen = nE-1;
N = nE;

phi1 =  evinfo.beta;
phi2 =  evinfo.alphas(1); 
eprior = evinfo.eprior; %[1, 100] 
logevidences = zeros(cycles,1);
newalphas = zeros(cycles,1);
newbetas = zeros(cycles,1);
errors = zeros(cycles,1);
bestbeta = phi1;
bestalpha = phi2;
bestlogevidence = 0;
bestres.bestant = 0;
bestres.besttour = [];
bestres.ETA = [];  
bestres.ants = []; 
% Lbest = evinfo.Lbest;
Lbest = inf;

% Note: |A| = SUMi( ln(Lambdi + alpha) ). Taken from part of equation 9.34

for k = 1:cycles

   
    anthessian = anthessianmatrix(ants, D, X, Y, ETA, nE);
    [~, evl] = eig(anthessian);
    evl = real(evl);                            % get the real number part of  the result.
    evl = evl .* (evl > 0);    
    
    mu = diag(evl);                             % expressed just before equation (9.38).Lamda(i) = Beta & mu(i) 
    L = phi2 * mu;
    
    gamma =  sum( L ./ (L + phi1));
    phi1 = 0.5 * gamma / eprior;                 % equation (9.41)
 
    logevidence1 = 0.5 * weightlen * log(phi1);  % portion of equation (9.33)
    
    logevidence2 = 0;    
    if( phi2 > 50.0)
        phi2 =  0.5 * ( N - gamma) / Error;       % equation (9.42)
           
        logevidence2 = 0.5 * N * log(phi2) - 0.5 * N * log(2 * pi);                         
        L = phi2 * mu;                    % ln|A|. taken from part of equation (9.34)       
    end
    
    sumlog =  L + phi1;
    sumlog =  (eps * (sumlog <= 0)) +  (sumlog .* (sumlog > 0)); % This sum is done this way to prevent 
                                                                 % negative or zero logs calculated in
                                                                 % the next line.    
    logevidence3 = 0.5 * sum( log(sumlog));
    
    sqerror = phi2 * Error + phi1 * eprior; 
    logevidence =  -sqerror - logevidence3 + logevidence1 + logevidence2;
                 
    paraminfo.alpha = phi2; % alpha;
    paraminfo.beta = phi1;  % beta;
    newalpha = phi2; 
    newbeta =  phi1;
      
    if evinfo.acomodel == 1 
       fprintf('\nasacotsp EF using new alpha and beta: %2.5f %2.5f\n', newalpha, newbeta); 
       res = asacotsp(paraminfo);
    elseif evinfo.acomodel == 2
       fprintf('\nacsacotsp EF using new alpha and beta: %2.5f %2.5f\n', newalpha, newbeta);  
       res = acsacotsp(paraminfo); 
    elseif evinfo.acomodel == 3
       fprintf('\nmmasacotsp EF using new alpha and beta: %2.5f %2.5f\n', newalpha, newbeta);  
       res = mmasacotsp(paraminfo); 
    elseif evinfo.acomodel == 4
       fprintf('\nrasacotsp EF using new alpha and beta: %2.5f %2.5f\n', newalpha, newbeta);  
       res = rasacotsp(paraminfo);  
    elseif evinfo.acomodel == 5
       fprintf('\nbwasacotsp EF using new alpha and beta: %2.5f %2.5f\n', newalpha, newbeta);  
       res = bwasacotsp(paraminfo);  
    elseif evinfo.acomodel == 6
       fprintf('\neasacotsp EF using new alpha and beta: %2.5f %2.5f\n', newalpha, newbeta);  
       res = easacotsp(paraminfo); 
    elseif evinfo.acomodel == 7
       fprintf('\nesacotsp ''17 EF using new alpha and beta: %2.5f %2.5f\n', newalpha, newbeta);  
       res = esacotsp(paraminfo); 
    elseif evinfo.acomodel == 8
       fprintf('\naacotsp ''22 EF using new alpha and beta: %2.5f %2.5f\n', newalpha, newbeta);  
       res = aaconctsp(paraminfo);    
    end  
    
    if evinfo.mode == 1
        
        if Lbest > res.besttour 
            if (res.besttour - benchmark) < 0                
               benchmark = res.besttour;
            else
                Error = res.besttour - benchmark;
            end
            
            if evinfo.acomodel ~= 7
               paraminfo.TAU = res.TAU;
               paraminfo.TOTAL = res.TOTAL;
            end 
            Lbest = res.besttour;           
            ants = res.ants;
            bestbeta = phi1; 
            bestalpha = phi2;
   
        else
           phi1 = bestbeta;
           phi2 = bestalpha;      
        end 
   
    elseif evinfo.mode == 2
        
        if Lbest > res.besttour 
           Lbest = res.besttour;
           Error = abs(res.besttour - benchmark);
           bestbeta = phi1; 
           bestalpha = phi2;
           bestlogevidence = logevidence;
           bestres = res;
           ants = res.ants;
           if evinfo.acomodel ~= 7
              paraminfo.TAU = res.TAU;
              paraminfo.TOTAL = res.TOTAL;
           end 
           
        else   
           phi1 = bestbeta;
           phi2 = bestalpha;
           logevidence = bestlogevidence; 
           res = bestres;
        end        
    else  
        
       Lbest = res.besttour;
       Error = abs(res.besttour - benchmark);
       if evinfo.acomodel ~= 7
          paraminfo.TAU = res.TAU;
          paraminfo.TOTAL = res.TOTAL;
       end       
    end
    
    logevidences(k) = logevidence;    
    newalphas(k) = phi2;
    newbetas(k) = phi1;
    errors(k) = Lbest;    
    respercycle(k) = res;                 %#ok<AGROW>
end 

end

function [hessian] = anthessianmatrix(ants, D, X, Y, ETA, nE)
   % The hessian matrix appears to indicate that the more ants
   % that traverse or move on edge e1, that's where the
   % minimum is. it does not nesserary translate to the smallest distance
   % in the graph or edge with the smallest distance.
   % there is other factors, like the greater the distance the smaller
   % the minimum.
  
   u = [1;1];
      
   hessian = zeros(nE);
   for k = 1:size(ants,1)         
       for p = 1:size(ants,2)-1    
    
           a = ants(k, p);
           b = ants(k, p + 1);
           if (a > b)
               b = a;
               a = ants(k, p + 1);
           end    
            
           dij = ETA(a,b);
           vi = [X(a); Y(a)];
           vj = [X(b); Y(b)];
           z = vi - vj;
           
           w = -2*u + ((2*dij) / D(a,b)) * ( u - (vi.*z)/D(a,b) ) + ...
                   (2*dij*(vj.*z))/ D(a,b);
           
           term = sqrt( w(1)^2 + w(2)^2);            
           if ~isinf(term) && ~isnan(term)
               hessian(a,b) = hessian(a,b) + term;
               hessian(b,a) = hessian(b,a) + term; 
           end      
        end
   end   
end
