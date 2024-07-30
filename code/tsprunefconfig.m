%% ------------------------------------------------------------------------
% tsprunefconfig.m is responsible for loading the data, initialising 
% the ants and other program parameters for all aco algirthms and  
% running them via the evidence framework.
%
% Use this to run various configurations of the same model.
%--------------------------------------------------------------------------

clc;
clear variables;

fprintf('---------------- Traveling Salesman Problem (TSP) --------------------\n');
fprintf('The program will be  loading data from the .mat files,\n');
fprintf('and then initialising the parameters and ants, ant size etc \n');
fprintf('and then invoke the relevant optimisation algorithm\n');
fprintf('-----------------------------------------------------------------------\n');

fprintf('Press any key to continue...\n\n');
pause;

directory  = '<directory of the .mat or tsp file>';
matfiles = dir(sprintf('%s/*.mat',directory));

%% ------------------------------------------------------------------------
%  RUN-TESTS for various ePRIORS for various datasets. 
%    
%-------------------------------------------------------------------------- 
acomodel = 1;
localsearch = 2;
filenumber = 1;

% -------------------------- LARGE TSPs CONFIGS --------------------------- 
%  Local search - 2 EPRIORS
%  OTHER-[40;80;55;40;70;50;50;50]; 
%  ESAC-[40;264;48;82;108;86;102;104];
%  AAC-[40;80;24;42;70;52;50;50];

%  Local search - 3 EPRIORS
%  OTHER-[40;80;16;34;70;50;50;50]; 
%  AS - [40;80;16;34;70;50;50;50]*8;
%  AAC-[40;80;24;42;70;52;50;50];
%  RAS-[48;80;8;25;128;18;18;34]*6;

EPRIORS = [40;80;55;40;70;50;50;50];

mode = 2;
evicycles = 2;     % OTHER-  evicycles:1 , ESACO - evicycles:2 |  3opt  OTHER - evicycles:2. RAS - evicycles:1
defaultalpha = 10; 
iteroffset = 0.25;  % 0.8; - AAC- 0.25;
modeltime = 0;
tspsize = 'MEDIUM';  % TSP size. LARGE or MEDIUM

%%-------------------------------------------------------------------------
% Fetch the samples generated by iRace and the results it calculated
%
%--------------------------------------------------------------------------
if strcmp(tspsize,'MEDIUM')
  iracedirectory =  '<directory of irace directory>/small-medium';
elseif strcmp(tspsize,'LARGE')    
  iracedirectory =  '<directory of irace directory>/large';    
end    

switch acomodel         
  case 1 % asacotspef 
        modeldirectory = sprintf('/AS/assamples%d%d.txt',localsearch,filenumber);
        arrayvecmodeldirectory = sprintf('/AS/asarrayvec%d%d.txt',localsearch,filenumber);
  case 2 % acsacotspef 
        modeldirectory = sprintf('/ACS/acssamples%d%d.txt',localsearch,filenumber);
        arrayvecmodeldirectory = sprintf('/ACS/acsarrayvec%d%d.txt',localsearch,filenumber);
  case 3 % mmasacotspef
        modeldirectory = sprintf('/MMAS/mmassamples%d%d.txt',localsearch,filenumber);
        arrayvecmodeldirectory = sprintf('/MMAS/mmasarrayvec%d%d.txt',localsearch,filenumber);
  case 4 % rasacotspef 
        modeldirectory = sprintf('/RAS/rassamples%d%d.txt',localsearch,filenumber);
        arrayvecmodeldirectory = sprintf('/RAS/rasarrayvec%d%d.txt',localsearch,filenumber);  
  case 6 % easacotspef
        modeldirectory = sprintf('/EAS/eassamples%d%d.txt',localsearch,filenumber);
        arrayvecmodeldirectory = sprintf('/EAS/easarrayvec%d%d.txt',localsearch,filenumber);
  case 7 % esacacotspef
        modeldirectory = sprintf('/ESAC/esacsamples%d%d.txt',localsearch,filenumber);     
        arrayvecmodeldirectory = sprintf('/ESAC/esacarrayvec%d%d.txt',localsearch,filenumber);
  case 8 % aacacotspef
        modeldirectory = sprintf('/AAC/aacsamples%d%d.txt',localsearch,filenumber);     
        arrayvecmodeldirectory = sprintf('/AAC/aacarrayvec%d%d.txt',localsearch,filenumber);      
  otherwise
       fprintf('\nIncorrect model select.\n');
end   

%% ------------------------------------------------------------------------
%  RUN-TESTS against various variables
%  rho = (0, 1]
%  antpop = [15, 80]; 
%  beta = random nos (1,10) randomnumbers(1,10, 10, 1) 
%  nnsize =  [20 30]
%  ranks = [6 100]
%  elitistants = [5 750]
%  cls = [2 50]
%  cand1 = [2 50]
%  cand2 = [2 10]
%% ------------------------------------------------------------------------

TOTALTIMES  = zeros(length(matfiles),1);

fileid  = fopen(sprintf('%s%s',iracedirectory, modeldirectory));
readline = fgetl(fileid);

maxsamplesize = 301;
nooffields = 11;     % rho, antpop, beta, nnsize, q0, ranks, elitistants cls cand1 cand2,alpha etc.
SAMPLES = zeros(maxsamplesize, nooffields); 
IRACEESULTS = zeros(maxsamplesize, 1); 
INSTANCES = strings(maxsamplesize, 1);
ANTSARRAYVEC = strings(maxsamplesize, 1);

samplesize = 0;

if acomodel == 8
   
   while ischar(readline)
       
      samplesize = samplesize + 1;
      sampleline = strtrim(readline);
      tokens = split(sampleline, ',');
      alphatoken = split(tokens{1},':');
      betatoken = split(tokens{2},':');
      rhotoken = split(tokens{3},':');
      antstoken = split(tokens{4},':');
      q0token = split(tokens{5},':');
      rankstoken = split(tokens{6},':');
      elitistantstoken = split(tokens{7},':');
      nnantstoken = split(tokens{8},':');     
      nsize = split(tokens{9},':');
      nsec = split(tokens{10},':');
      resultstoken = split(tokens{13},':');
      instancetoken = split(tokens{14},':');     
      
      SAMPLES(samplesize, 1) = str2double(strtrim(rhotoken{2}));
      SAMPLES(samplesize, 2) = str2double(strtrim(antstoken{2}));
      SAMPLES(samplesize, 3) = str2double(strtrim(betatoken{2}));
      SAMPLES(samplesize, 4) = str2double(strtrim(nnantstoken{2}));
      SAMPLES(samplesize, 5) = str2double(strtrim(q0token{2}));
      SAMPLES(samplesize, 6) = str2double(strtrim(rankstoken{2}));
      SAMPLES(samplesize, 7) = str2double(strtrim(elitistantstoken{2}));     
      SAMPLES(samplesize, 8) = 0;
      SAMPLES(samplesize, 9) = str2double(strtrim(nsize{2}));
      SAMPLES(samplesize, 10) = str2double(strtrim(nsec{2}));     
      SAMPLES(samplesize, 11) = str2double(strtrim(alphatoken{2}));
     
      IRACEESULTS(samplesize, 1) = str2double(strtrim(resultstoken{2}));
      INSTANCES(samplesize, 1) = strtrim(instancetoken{2});     
      readline = fgetl(fileid);
   end 
    
else
    
   while ischar(readline)
       
      samplesize = samplesize + 1;
      sampleline = strtrim(readline);
      tokens = split(sampleline, ',');
      alphatoken = split(tokens{1},':');
      betatoken = split(tokens{2},':');
      rhotoken = split(tokens{3},':');
      antstoken = split(tokens{4},':');
      q0token = split(tokens{5},':');
      rankstoken = split(tokens{6},':');
      elitistantstoken = split(tokens{7},':');
      nnantstoken = split(tokens{8},':');     
      cltoken = split(tokens{9},':');
      cand1token = split(tokens{10},':');
      cand2token = split(tokens{11},':');
      resultstoken = split(tokens{14},':');
      instancetoken = split(tokens{15},':');     
      SAMPLES(samplesize, 1) = str2double(strtrim(rhotoken{2}));
      SAMPLES(samplesize, 2) = str2double(strtrim(antstoken{2}));
      SAMPLES(samplesize, 3) = str2double(strtrim(betatoken{2}));
      SAMPLES(samplesize, 4) = str2double(strtrim(nnantstoken{2}));
      SAMPLES(samplesize, 5) = str2double(strtrim(q0token{2}));
      SAMPLES(samplesize, 6) = str2double(strtrim(rankstoken{2}));
      SAMPLES(samplesize, 7) = str2double(strtrim(elitistantstoken{2}));     
      SAMPLES(samplesize, 8) = str2double(strtrim(cltoken{2}));
      SAMPLES(samplesize, 9) = str2double(strtrim(cand1token{2}));
      SAMPLES(samplesize, 10) = str2double(strtrim(cand2token{2}));     
      SAMPLES(samplesize, 11) = str2double(strtrim(alphatoken{2}));
     
      IRACEESULTS(samplesize, 1) = str2double(strtrim(resultstoken{2}));
      INSTANCES(samplesize, 1) = strtrim(instancetoken{2});     
      readline = fgetl(fileid);
   end

end
fclose(fileid);

fileid  = fopen(sprintf('%s%s',iracedirectory, arrayvecmodeldirectory),'r');
readline = fgetl(fileid);
vecsize = 0;
while ischar(readline)
    
     vecsize = vecsize + 1;
     sampleline = strtrim(readline);
     ANTSARRAYVEC(vecsize, 1) = sampleline;     
     readline = fgetl(fileid);
end    
fclose(fileid);

MINBETAS = zeros(samplesize, 1);
MINERRORS = zeros(samplesize, 1);
FILESAMPLEEVIDENCEERROR = zeros(samplesize, 3);                                          

for in = 1 : samplesize
        
    for fileindex = 1 : length(matfiles)       
        if matfiles(fileindex).name == INSTANCES(in)
           filename  =  matfiles(fileindex).name;
           load(sprintf('%s\\%s',matfiles(fileindex).folder, matfiles(fileindex).name));
           break;
        end        
    end    
    
    paraminfo.X = X;
    paraminfo.Y = Y;
    paraminfo.nE = nE;
    paraminfo.Distance = Distance; 
    paraminfo.VecDistance = reshape(Distance,[],1);
    paraminfo.Q = 1;
    paraminfo.maxiteration = maxiteration * iteroffset; 
    paraminfo.localsearch = localsearch;
    paraminfo.q0 = 0;
    paraminfo.ranks = 0;
    
    paraminfo.alpha = SAMPLES(in, 11); 
    if acomodel == 7
       if paraminfo.alpha == 0
           paraminfo.alpha = 0.9;
       end    
  
       paraminfo.maxclcount = SAMPLES(in, 8);
       paraminfo.numberofcandidatesconst = SAMPLES(in, 9);
       paraminfo.numberofcandidates  = SAMPLES(in, 10);
    elseif acomodel == 8  
       paraminfo.nsize = SAMPLES(in, 9);
       paraminfo.nsectors = SAMPLES(in, 10);
       
    else  
       if paraminfo.alpha == 0
           paraminfo.alpha = eps;
       end        
    end    
    
    paraminfo.rho = SAMPLES(in, 1); 
    paraminfo.antpop = SAMPLES(in, 2);
    paraminfo.beta = SAMPLES(in, 3); 
    paraminfo.nnsize = SAMPLES(in, 4);
    
    paraminfo.q0 = SAMPLES(in ,5);
    paraminfo.ranks = SAMPLES(in ,6);
    paraminfo.elitistants = SAMPLES(in, 7);
    paraminfo.fullpathname = fullpathname;
    paraminfo.nearneighbours = nnlist(Distance, nE, paraminfo.nnsize); 
          
    evinfo.cycles = evicycles; 
    evinfo.alphas = paraminfo.alpha;
    evinfo.beta = paraminfo.beta;   
    evinfo.X = X ;
    evinfo.Y = Y;
    evinfo.Distance = Distance;         
    evinfo.nE = nE;       
    evinfo.acomodel = acomodel; 
    evinfo.eprior = EPRIORS(fileindex);  
    evinfo.mode = mode;

    strtokens = split(ANTSARRAYVEC(in, 1)," ");     
    doublevec = str2double(strtokens)'; 
    evinfo.ants =  reshape(doublevec,paraminfo.antpop, nE); % res.ants;
    evinfo.Lbest = IRACEESULTS(in, 1); 
    evinfo.benchmark = benchmarks(fileindex);
    evinfo.ETA = 1./(Distance + 0.1); 
    evinfo.Error = abs(evinfo.Lbest - benchmarks(fileindex)) / 2;

    tic;
    fprintf('\nBest tour before the rho evidence: %d. (Sample %d of %d), instance: %s\n', evinfo.Lbest, in, samplesize, filename);    
    [newalphas, newbetas, logevidences, errors, respercycle] = evidenceacotsp(evinfo, paraminfo);
    [minerror, minindex] = min(errors);
    TOTALTIMES(fileindex) =  TOTALTIMES(fileindex) + toc;
    
    MINBETAS(in) = newbetas(minindex);
    MINERRORS(in) = minerror; 
    
    FILESAMPLEEVIDENCEERROR(in, 1) = fileindex;
    FILESAMPLEEVIDENCEERROR(in, 2) = evinfo.Lbest;
    FILESAMPLEEVIDENCEERROR(in, 3) = minerror; 
       
end    

resultsdirectiory = '<directory of results>';

switch acomodel         
  case 1
       % ================================  asacotspef ====================
       savefile = sprintf('%sASACOTSP%d',tspsize,filenumber); 
  case 2  
       % ================================  acsacotspef ===================
       savefile = sprintf('%sACSACOTS%d',tspsize,filenumber); 
  case 3
       % ================================  mmasacotspef ==================
       savefile = sprintf('%sMMASACOTSP%d',tspsize,filenumber);  
  case 4
       % ================================  rasacotspef ===================
       savefile = sprintf('%sRASACOTSP%d',tspsize,filenumber);  
   case 6 
       % ================================  easacotspef ===================
       savefile = sprintf('%sEASACOTSP%d',tspsize,filenumber); 
  case 7 
       % ================================  esacotspef ===================
       savefile = sprintf('%sESAACO17TSP%d',tspsize,filenumber);      
  case 8 
       % ================================  aaconctsp ===================
       savefile = sprintf('%sAACO22TSP%d',tspsize,filenumber);           
  otherwise
       fprintf('\nIncorrect model select.\n');
end   

save(sprintf('%s/%s.mat',resultsdirectiory, savefile), 'nE','maxiteration','antpop','acomodel','Distance', 'X','Y','type','weighttype','matfiles',...
                                                      'evicycles','name','modeltime','samplesize','localsearch','TOTALTIMES',...
                                                      'MINBETAS','MINERRORS','FILESAMPLEEVIDENCEERROR', '-v7.3');
                                                                 
                                                               