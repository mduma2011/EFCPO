%% -------------------------------------------------------------
% tspun.m is responsible for loading the data, initialising the
% ants and other program parameters   
%
%---------------------------------------------------------------

clc;
clear variables;

fprintf('---------------- Traveling Salesman Problem (TSP) --------------------\n');
fprintf('The program will be  loading data from the .mat files,\n');
fprintf('and then initialising the parameters and ants \n');
fprintf('and then invoke the relevant optimisation algorithm\n');
fprintf('-----------------------------------------------------------------------\n');

fprintf('Press any key to continue...\n\n');
pause;

directory  = '<directory of the .mat or .tsp file>';
matfiles = dir(sprintf('%s/*.mat',directory));
RES = zeros(8, 8);
TIMES = zeros(8, 8);

count = 1;
for fileindex = 1: length(matfiles)
    
    filename  =  matfiles(fileindex).name;
    load(sprintf('%s\\%s',matfiles(fileindex).folder, matfiles(fileindex).name));
    
    paraminfo.X = X;
    paraminfo.Y = Y;
    paraminfo.nE = nE;
    paraminfo.Distance = Distance;   
    paraminfo.VecDistance = reshape(Distance,[],1);
    paraminfo.maxiteration = maxiteration; % maxiteration*20;
    paraminfo.antpop = antpop;
    paraminfo.Q = 1;
    paraminfo.localsearch = 2;  % 0- no search: 2: 2-opt ,3: 3-opt
    
    % ================================  asacotsp ==========================
    
    paraminfo.alpha = 1;
    paraminfo.beta = 7.175416; % 7.431189;  
    paraminfo.rho = 0.831117; % 0.4647;
    paraminfo.nnsize = 54;
    paraminfo.antpop = 100;
    
    if paraminfo.nnsize >= nE
       paraminfo.nnsize = 20; 
    end    
    
    paraminfo.nearneighbours = nnlist(Distance, nE, paraminfo.nnsize); 
    
    fprintf('\n1. Executing the TSP asacotspef against file: %s\n', filename);

    tic;
    res1 = asacotsp(paraminfo);
    elapsedTime1 = toc; 
         
    % ================================  acsacotsp =========================
    
    paraminfo.alpha = 2.67;
    paraminfo.beta = 7.16;  
    paraminfo.rho = 0.5;
    paraminfo.nnsize = 49;
    paraminfo.q0 = 0.05;
    paraminfo.nearneighbours = nnlist(Distance, nE, paraminfo.nnsize); 
      
    fprintf('\n2. Executing the TSP acsacotspef against file: %s\n', filename);
    
    tic;
    res2 = acsacotsp(paraminfo);
    elapsedTime2 = toc;
         
    % ================================  mmasacotsp =========================
    
    paraminfo.alpha = 4.08;
    paraminfo.beta = 5.58; 
    paraminfo.rho = 0.47;
    paraminfo.nnsize = 73;
     
    if paraminfo.nnsize >= nE
       paraminfo.nnsize = 20; 
    end    
    paraminfo.q0 = 0;
    paraminfo.nearneighbours = nnlist(Distance, nE, paraminfo.nnsize); 
     
    fprintf('\n3. Executing the TSP mmasacotsp against file: %s\n', filename);
    
    tic;
    res3 = mmasacotsp(paraminfo);
    elapsedTime3 = toc;     
     
    % ================================  rasacotsp =========================
    
    paraminfo.alpha = 3.85;
    paraminfo.beta =  5.32; 
    paraminfo.rho = 0.76;
    paraminfo.nnsize = 31;
    paraminfo.q0 = 0;
    paraminfo.ranks = 71;
    paraminfo.nearneighbours = nnlist(Distance, nE, paraminfo.nnsize); 
     
    fprintf('\n4. Executing the TSP rasacotspef against file: %s\n', filename);
    
    tic;
    res4 = rasacotsp(paraminfo);
    elapsedTime4 = toc; 
    
            
    % ================================  easacotsp =========================
    
    paraminfo.alpha = 3.85;
    paraminfo.beta = 5.32; 
    paraminfo.rho = 0.76;
    paraminfo.nnsize = 31;
     
    if paraminfo.nnsize >= nE
       paraminfo.nnsize = 20; 
    end    
    paraminfo.q0 = 0;
    paraminfo.ranks = 0;
    paraminfo.elitistants = 51;
    paraminfo.nearneighbours = nnlist(Distance, nE, paraminfo.nnsize); 
     
    fprintf('\n6. Executing the TSP rasacotspef against file: %s\n', filename);
    
    tic;
    res6 = easacotsp(paraminfo);
    elapsedTime6 = toc; 

    % ================================  esacotsp =========================

    paraminfo.alpha = 3.85;
    paraminfo.beta = 5.32; 
    paraminfo.rho = 0.76;
    paraminfo.nnsize = 31;
    paraminfo.maxclcount = 10;
    paraminfo.numberofcandidatesconst = 20;
    paraminfo.numberofcandidates = 4;
    paraminfo.fullpathname = fullpathname;  
       
    if paraminfo.nnsize >= nE
        paraminfo.nnsize = 20; 
    end    
       
    fprintf('\n7. Executing the TSP esacotsp ''17 against file: %s\n', filename);
       
    tic;
    res7 = esacotsp(paraminfo);
    elapsedTime7 = toc;
    
    % ================================ aaconctsp =========================
    
     paraminfo.alpha = 1.53;
     paraminfo.beta = 1.281; 
     paraminfo.rho = 0.1;  
     paraminfo.nnsize = 20;
     paraminfo.nsize = 8;
     paraminfo.nsectors = 3;
    
     if paraminfo.nnsize >= nE
        paraminfo.nnsize = 20; 
     end    
    
     paraminfo.nearneighbours = nnlist(Distance, nE, paraminfo.nnsize); 
    
     fprintf('\n1. Executing the TSP aaconctsp against file: %s\n', filename);
     
     tic;
     res8 = aaconctsp(paraminfo);
     elapsedTime8 = toc; 
     
    % ================================== end ==============================
        
     RES(1, count) = res1.besttour; 
     RES(2, count) = res2.besttour;
     RES(3, count) = res3.besttour;
     RES(4, count) = res4.besttour;
     RES(5, count) = res5.besttour;
     RES(6, count) = res6.besttour;
     RES(7, count) = res7.besttour;
     RES(8, count) = res8.besttour;
     
     TIMES(1, count) = elapsedTime1;
     TIMES(2, count) = elapsedTime2;
     TIMES(3, count) = elapsedTime3;
     TIMES(4, count) = elapsedTime4;
     TIMES(5, count) = elapsedTime5;
     TIMES(6, count) = elapsedTime6;
     TIMES(7, count) = elapsedTime7;
     TIMES(8, count) = elapsedTime8;
     
     count = count + 1;
end    
fprintf('Finished. Press any key to end...\n\n');
pause;