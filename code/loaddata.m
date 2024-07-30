%%-------------------------------------------------------------------------
%   Ant colony Optimisation algorithm based on the  
%   implementation by T. Stuetzle 
%
%  
%    
%%-------------------------------------------------------------------------

%% ------------------------------------------------------------------------
%  1. Read the data from the files
% -------------------------------------------------------------------------
clear;
clc;
fprintf('---------------- Traveling Salesman Problem (TSP) --------------------\n');
fprintf('The program will be reading or loading data from text or tsp files,\n');
fprintf('and then populating the data into nE variable and Distance matrices,\n');
fprintf('and thereafter, saving them into the relevant .mat files.\n');
fprintf('-----------------------------------------------------------------------\n');
fprintf('Press any key to continue...\n');
pause;

 
directory  = '<directory of the tsp file>';
          
filenames =  ["pr1002.tsp","pr2392.tsp","rat783.tsp","rl1323.tsp",...
              "rl1889.tsp","u1432.tsp","u1817.tsp","u2152.tsp",...
              "pcb3038.tsp","rl5934.tsp","usa13509.tsp"];

benchmarks = [259045,378032,8806,270199,...
              316536,152970,57201,64253,...
              137694,554070,19947008];

fileindices = [1, 2, 3, 4, 5, 6, 7, 8,...
               9, 10, 11];                
                        
maxiterations = [100,100,100,100,100,100,100,100,100,...
                100,100];
            
antpops = [25,25,25,25,25,25,25,25,...
           25,25,25];           
       
                        

for index = 1:length(filenames)
    
 fullpathname = sprintf('%s/%s',directory, filenames(index));   
 fileid  = fopen(fullpathname);
 maxiteration = maxiterations(index);
 antpop = antpops(index);
      
 readline = fgetl(fileid);
 readcoordinates = false;
 node = 1;
 while ischar(readline)
       
     field = strtrim(readline);
     
     if startsWith(field, "NAME") || startsWith(field, "NAME:")
         strplit =  split(field,":");
         name = strtrim(strplit{2});
         name  = string(upper(name));
     
     elseif  startsWith(field, "COMMENT") || startsWith(field, "COMMENT:")
         strplit =  split(field,":");
         comment = strtrim(strplit{2});
         comment  = string(upper(comment));
         
     elseif  startsWith(field, "TYPE") || startsWith(field, "TYPE:")
         strplit =  split(field,":");
         type = strtrim(strplit{2});
         type  = string(upper(type));
     
     elseif  startsWith(field, "DIMENSION") || startsWith(field, "DIMENSION:")
         strplit =  split(field,":");
         dim = strtrim(strplit{2});
         dim  = string(upper(dim));
         nE =  str2double(dim);
         X = zeros(nE, 1);
         Y = zeros(nE, 1);     
         
      elseif  startsWith(field, "EDGE_WEIGHT_TYPE") || startsWith(field, "EDGE_WEIGHT_TYPE:")
         strplit =  split(field,":");
         weighttype = strtrim(strplit{2});
         weighttype  = string(upper(weighttype));   
     end
     
     if strcmp(field,'EOF') || strcmp(field,'TOUR_SECTION')
        readcoordinates = false; 
     end    
         
     if readcoordinates           
        splits = strsplit(field);
        x = str2double(splits{1,2}); 
        y = str2double(splits{1,3}); 
            
        X(node) = x;
        Y(node) = y;
        node = node + 1;
     end
    
     if strcmp(field,'NODE_COORD_SECTION') ||...
        strcmp(field,'DISPLAY_DATA_SECTION')
        readcoordinates = true;                   
     end  
     readline = fgetl(fileid);
  end         
 
 fclose(fileid);
 
 switch(weighttype)
  
     case 'EUC_2D'
          Distance =  euc_distance(X, Y, nE);
     case 'CEIL_2D'
          Distance =  ceil_distance(X, Y, nE);
     case 'GEO'
          Distance =  geo_distance(X, Y, nE);
     case 'ATT'
          Distance =  att_distance(X, Y, nE);
     otherwise
        fprintf('EDGE_WEIGHT_TYPE: %s was not implemented\n\n',weighttype);
 end       
 
  % 2. Save the data the data to a .mat file.
  %    for now we will save the number of edges (nE) and the Distance matrix
  savefile = name;
  save(sprintf('%s/%s.mat',directory, savefile), 'nE','maxiteration','antpop','Distance', 'X','Y','name','comment',...
                                                 'weighttype','type','benchmarks','fullpathname','-v7.3'); 
end

fprintf('\n---------- DONE loading data for the Traveling Salesman Problem -----------\n');
fprintf('Press any key to continue...\n');
pause;

function Distance = euc_distance(X, Y, nE)
% EUC_2D Function
% Computes the euclidean distance between two nodes
  Distance = zeros(nE, nE);  
  for i = 1:nE -1
     for j = i + 1:nE
         Distance(i,j) = sqrt((X(i)-X(j))^2 + (Y(i)-Y(j))^2) + 0.5;
         Distance(j,i) = Distance(i,j);
     end   
  end
  Distance = floor(Distance);
end

function Distance = ceil_distance(X, Y, nE)
% CEIL_2D Function
% Computes the ceiling distance between two nodes rounded to next  integer
  Distance = zeros(nE, nE);  

  for i = 1:nE -1
     for j = i + 1:nE
         Distance(i,j) = sqrt((X(i)-X(j))^2 + (Y(i)-Y(j))^2);
         Distance(j,i) = Distance(i,j);
     end   
  end
  Distance = ceil(Distance);
end

function Distance = geo_distance(X, Y, nE)
% GEO Function
% Compute geometric distance between two nodes rounded to next integer
  Distance = zeros(nE, nE);
  
  for i = 1:nE -1
     for j = i + 1:nE
         
         x1 = X(i);  y1 = Y(i);
         x2 = X(j);  y2 = Y(j);
         
         deg = floor(x1);
         min  = x1 - deg;
         lati = pi * (deg + 0.5 * min / 3.0) / 180.0;         
         deg = floor(x2);
         min  = x2 - deg;
         latj = pi * (deg + 0.5 * min / 3.0) / 180.0;
         
         deg = floor(y1);
         min  = y1 - deg;
         longi = pi * (deg + 0.5 * min / 3.0) / 180.0;         
         deg = floor(y2);
         min  = y2 - deg;
         longj = pi * (deg + 0.5 * min / 3.0) / 180.0;
         
         c1 = cos (longi - longj);
         c2 = cos (lati - latj);
         c3 = cos (lati + latj);
          
         Distance(i,j) = floor( (6378.388 * acos (0.5 * ((1.0 + c1) * c2 - (1.0 - c1) * c3)) + 1.0));
         Distance(j,i) = Distance(i,j);
     end   
  end
  

end

function Distance = att_distance(X, Y, nE)
% ATT Function
% Compute ATT distance between two nodes rounded to next integer
  Distance = zeros(nE, nE);  
  
  for i = 1:nE -1
     for j = i + 1:nE
        
         Rij = sqrt( ((X(i)-X(j))^2 + (Y(i)-Y(j))^2)/10.0 );
         Tij  =  floor(Rij);
         if Tij < Rij
            Distance(i,j) = Tij + 1;  
         else
            Distance(i,j) = Tij; 
         end    
         Distance(j,i) = Distance(i,j);
     end   
  end

end