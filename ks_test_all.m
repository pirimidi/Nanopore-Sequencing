%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: February 4, 2016.
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: Given 4 'stats.mat' structure of static capture vs. background
% data sets, this program iterates through all of them and generates the
% p-value for each pair (A,C,T,G) using the 2D Kolmogorov-Smirnov test for:
%
% (1) DE := dwell time of all filtered events (exp)
% (2) DB := dwell time of all filtered events (bgr)
% (3) DEB := dwell time of events in tag cature band (exp)
% (4) DBB := dwell time of events in tag capture band (bgr)
% (5) DET := dwell time of events in tag cature zone (exp)
% (6) DBT := dwell time of events in tag capture zone (bgr)
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function ks_test_all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         EVENT ANALYSIS STARTUP                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all')

fprintf('\n');
disp('--> KS test for static capture versus background start');
fprintf('\n');

% Set default number formatting.
format short;

% Define current working directory.
work_dir = pwd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         DATA RETRIEVAL SECTION                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> DATA RETRIVAL SECTION');

% Navigate to 'stats' data directory.
if ~exist('stats', 'dir')
  mkdir('stats');
end

cd 'stats';

% Read in all 'stats' text file names one-by-one.
list = dir('stats_*');

% Define array container for dwell time and normalized current of all 
% filtered events (exp).
S = struct('DE', {}, 'CE', {}, 'DB', {}, 'CB', {}, ...
           'DEB', {}, 'CEB', {}, 'DBB', {}, 'CBB', {}, ...
           'DET', {}, 'CET', {}, 'DBT', {}, 'CBT', {});

for i = 1:length(list)

  disp(['--> Processing structure: ', list(i).name]);  

  % Load the stats data structure into workspace.
  load(list(i).name);

  % Generate container arrays for all filtered events (exp).
  S(i).DE = dwellt_all_exp;	
  S(i).CE = ncurr_all_exp;
  
  % Generate container arrays for all filtered events (bgr).
  S(i).DB = dwellt_all_bgr;
  S(i).CB = ncurr_all_bgr;
  
  % Generate container arrays for events in tag capture band(exp).
  S(i).DEB = dwellt_tcb_exp;
  S(i).CEB = ncurr_tcb_exp;  
  
  % Generate container arrays for events in tag capture band (bgr).
  S(i).DBB = dwellt_tcb_bgr;
  S(i).CBB = ncurr_tcb_bgr;    

  % Generate container arrays for events in tag capture zone (exp).
  S(i).DET = dwellt_tcz_exp;
  S(i).CET = ncurr_tcz_exp;  
  
  % Generate container arrays for events in tag capture zone (bgr).
  S(i).DBT = dwellt_tcz_bgr;
  S(i).CBT = ncurr_tcz_bgr;    

end

fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     CALCULATE P-VALUES USING KS TEST                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Dwell time for all data (DE/DB).
disp('--> 1. DWELL TIME P-VALUES USING KS TEST: ALL DATA (DE/DB)');

for i = 1:length(list)
    for j = 1:length(list)

        disp(['--> Matrix index: ', num2str(i), '-', num2str(j)]);
   
        % Compute p-value of KS test.
        [h1(i,j), p1(i,j), ks2stat1(i,j)] = kstest2(S(i).DE, S(j).DB);
        %[h, p] = kstest2(S(i).DE, S(j).DB)
                
    end
end

h1
p1
ks2stat1

% 2. Dwell time for tag capture zone only (DET/DBT).
disp('--> 2. DWELL TIME P-VALUES USING KS TEST: TAG CAPTURE ZONE ONLY (DET/DBT)');

for i = 1:length(list)
    for j = 1:length(list)

        disp(['--> Matrix index: ', num2str(i), '-', num2str(j)]);
   
        % Compute p-value of KS test.
        [h2(i,j), p2(i,j), ks2stat2(i,j)] = kstest2(S(i).DET, S(j).DBT);
        [h, p] = kstest2(S(i).DET, S(j).DBT)
                
    end
end

h2
p2
ks2stat2

% 3. Dwell time for tag capture band only (DEB/DBB).
disp('--> 3. DWELL TIME P-VALUES USING KS TEST: TAG CAPTURE BAND ONLY (DEB/DBB)');

for i = 1:length(list)
    for j = 1:length(list)

        disp(['--> Matrix index: ', num2str(i), '-', num2str(j)]);
   
        % Compute p-value of KS test.
        [h3(i,j), p3(i,j), ks2stat3(i,j)] = kstest2(S(i).DEB, S(j).DBB);
        [h, p] = kstest2(S(i).DEB, S(j).DBB)
                
    end
end

h3
p3
ks2stat3

% 4. Normalized current for all data (CE/CB).
disp('--> 4. CURRENT P-VALUES USING KS TEST: ALL DATA (CE/CB)');

for i = 1:length(list)
    for j = 1:length(list)

        disp(['--> Matrix index: ', num2str(i), '-', num2str(j)]);
   
        % Compute p-value of KS test.
        [h4(i,j), p4(i,j), ks2stat4(i,j)] = kstest2(S(i).CE, S(j).CB);
        %[h, p] = kstest2(S(i).CE, S(j).CB)
                
    end
end

h4
p4
ks2stat4

% 5. Normalized current for tag capture band only (CEB/CBB).
disp('--> 5. CURRENT P-VALUES USING KS TEST: TAG CAPTURE BAND ONLY (CEB/CBB)');

for i = 1:length(list)
    for j = 1:length(list)

        disp(['--> Matrix index: ', num2str(i), '-', num2str(j)]);
   
        % Compute p-value of KS test.
        [h5(i,j), p5(i,j), ks2stat5(i,j)] = kstest2(S(i).CEB, S(j).CBB);
        %[h, p] = kstest2(S(i).CET, S(j).CBT)
                
    end
end

h5
p5
ks2stat5

% 6. Normalized current for tag capture zone only (CET/CBT).
disp('--> 6. CURRENT P-VALUES USING KS TEST: TAG CAPTURE ZONE ONLY (CET/CBT)');

for i = 1:length(list)
    for j = 1:length(list)

        disp(['--> Matrix index: ', num2str(i), '-', num2str(j)]);
   
        % Compute p-value of KS test.
        [h6(i,j), p6(i,j), ks2stat6(i,j)] = kstest2(S(i).CET, S(j).CBT);
        %[h, p] = kstest2(S(i).CET, S(j).CBT)
                
    end
end

h6
p6
ks2stat6

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             KS TEST FINALIZATION                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> KS test for static capture versus background end');
fprintf('\n');
