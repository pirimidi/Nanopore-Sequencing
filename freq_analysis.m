%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: August 26, 2015.
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: Given a 'var.mat' structure of a 'vardat' run for N data sets,
% this program iterates through all of them and extracts filter level 2 
% (all clean pores) event information of current blockade level and 
% corresponding dwell time pairs, then generates a confusion matrix answering 
% the following question:
%
% "How many times are we in the signal band(s) for all capture and back-
% ground events?" (Answer is given as a % of total events.)
%
% Input arguments:
%
% (1) filter_e := filter level for pores in experimrntal data [1-4]
% (2) filter_c := filter level for pores in background data [1-4]
% (3) order_e := stack order, i.e., experimental data sets to access
% (4) order_c := stack order, i.e., background data sets to access
% (5) tag_capture := lower and upper tag capture limits, an array of {GACT}
% (6) dwell_cutoff := background dwell time cutoff (in sec), an array
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function [A] = freq_analysis(filter_e, filter_b, order_e, order_b, ...
                             tag_capture, dwell_cutoff)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         ERROR ANALYSIS STARTUP                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all')

fprintf('\n');
disp('--> Frequency analysis start');
fprintf('\n');

% Set default number formatting.
format short;

% Define current working directory.
work_dir = pwd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  'EXPERIMENT' ERROR ANALYSIS SECTION                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> EXPERIMENT EVENT STATISTICS SECTION');

% Navigate to 'experiment' data directory.
if ~exist('experiment', 'dir')
  mkdir('experiment');
end

cd 'experiment';

% Load the event data structure into workspace.
load('var.mat');

% Read in all 'experiment' statistic text files one-by-one.
list = dir('data*_vardat');

% Define array container for dwell time data.
DE = []; 

% Define array container for current blockade level data.
CE = [];

for j = 1:length(order_e)

  disp(['--> Processing file: ', list(order_e(j)).name]);  

  % Access data stack for dataset 'j'.
  stack = var.stack(order_e(j));	

  % Access data structure containing filter level 2 (all clean single 
  % channels) event information.
  cap = stack.stat(filter_e).cap;

  % Retrieve all current blockade level and corresponding dwell time pairs.
  dwell = cap.dwell(:, :);
  imed = cap.imed(:, :);

  % Filter out NaNs and diplay 2D scatter plot and histogram with a 
  % distribution fit.
  iD = ~isnan(dwell);
  de = dwell(iD(:));
  ce = imed(iD(:));

  % Obtain size of maximum events.
  [r, c] = size(de);

  % Display the number of maximum events captured by a single pore.
  disp(['--> Maximum number of events: ', num2str(r)]);
  fprintf('\n');

  % Collect all filtered events and update storage array.
  DE = [DE; de];
  CE = [CE; ce];

end

scatter(DE, CE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       FILTER OUT UPPER BAND                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> TAG CAPTURE EVENT STATISTICS SECTION');

% Define array container for dwell time and current blockade level data of 
% 'experiment' data set.
DEU = []; 
CEU = [];

% Counter of tag capture events.
q = 1;

% Filter out upper events in the 'tag capture' for 'experiment' data set.
for w = 1:length(DE)
    
    % If event is in not lower tag capture.
    if not(CE(w) < 0.32 && DE(w) > 0.01) 

            % Remember these tag capture events.
            DEU(q) = DE(w);
            CEU(q) = CE(w);

            % Increment tag capture event counter.
            q = q+1;
    end  
end

% Obtain size of all filtered events.
[r, c] = size(DEU);

% Display the number of maximum events captured by a single pore.
disp(['--> Number of all filtered events (exp): ', num2str(r)]);
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  'BACKGROUND' ERROR ANALYSIS SECTION                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> BACKGROUND EVENT STATISTICS SECTION');

% Navigate to working directory.
cd(work_dir);

% Navigate to 'background' data directory.
if ~exist('background', 'dir')
  mkdir('background');
end

cd 'background';

% Load the event data structure into workspace.
load('var.mat');

% Read in all 'background' statistic text files one-by-one.
list = dir('data*_vardat');

% Define array container for dwell time data.
DB = []; 

% Define array container for current blockade level data.
CB = [];

for l = 1:length(order_b)

  disp(['--> Processing file: ', list(order_b(l)).name]);  

  % Access data stack for dataset 'l'.
  stack = var.stack(order_b(l));	

  % Access data structure containing filter level 2 (all clean single 
  % channels) event information.
  cap = stack.stat(filter_b).cap;

  % Retrieve all current blockade level and corresponding dwell time pairs.
  dwell = cap.dwell(:, :);
  imed = cap.imed(:, :);

  % Filter out NaNs and diplay 2D scatter plot and histogram with a 
  % distribution fit.
  iD = ~isnan(dwell);
  db = dwell(iD(:));
  cb = imed(iD(:));

  % Obtain size of maximum events.
  [r, c] = size(db);

  % Display the number of maximum events captured by a single pore.
  disp(['--> Maximum number of events: ', num2str(r)]);
  fprintf('\n');

  % Collect all filtered events and update storage array.
  DB = [DB; db];
  CB = [CB; cb];

end

% Obtain size of all filtered events.
[r, c] = size(DB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  'TAG CAPTURE' ERROR ANALYSIS SECTION                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> TAG CAPTURE ERROR ANALYSIS SECTION - EXP');

% Define array container for dwell time and current blockade level data of 
% 'experiment' data sets = {G, A, C, T}.
DETG = []; CETG = [];
DETA = []; CETA = [];
DETC = []; CETC = [];
DETT = []; CETT = [];

% Counter of events in tag capture box.
u = 0;

% Counter of G/A/C/T-tag capture events.
ge = 1; ae = 1; cc = 1; te = 1;

% Filter out events in the 'tag capture' for 'experiment' data set.
for m = 1:length(DEU)

    % If event is in dwell time cutoff range (x) as well as in capture 
    % range (y), i.e., in tag capture box. 
    if (dwell_cutoff(1) < DEU(m) && DEU(m) < dwell_cutoff(2)) && ...
       (tag_capture(1) < CEU(m) && CEU(m) < tag_capture(8))  
   
        % Increment events in tag capture box.
        u = u+1;
    
        % If event is in G-tag capture range.
        if tag_capture(1) < CEU(m) && CEU(m) < tag_capture(2) 
            
            % Remember these G-tag capture events.
            DETG(ge) = DEU(m);
            CETG(ge) = CEU(m);
            
            % Increment boxed-in G-tag capture event counter.
            ge = ge+1;
            
        end
        
        % If event is in A-tag capture range.
        if tag_capture(3) < CEU(m) && CEU(m) < tag_capture(4) 
            
            % Remember these A-tag capture events.
            DETA(ae) = DEU(m);
            CETA(ae) = CEU(m);
            
            % Increment boxed-in A-tag capture event counter.
            ae = ae+1;
            
        end
        
        % If event is in C-tag capture range.
        if tag_capture(5) < CEU(m) && CEU(m) < tag_capture(6) 
            
            % Remember these C-tag capture events.
            DETC(cc) = DEU(m);
            CETC(cc) = CEU(m);
            
            % Increment boxed-in C-tag capture event counter.
            cc = cc+1;
            
        end
        
        % If event is in T-tag capture range.
        if tag_capture(7) < CEU(m) && CEU(m) < tag_capture(8) 
            
            % Remember these T-tag capture events.
            DETT(te) = DEU(m);
            CETT(te) = CEU(m);
            
            % Increment boxed-in T-tag capture event counter.
            te = te+1;
            
        end        
    end  
end

% Obtain size of all filtered and tag capture events.
[re, cw] = size(DE);
[retg, cetg] = size(DETG);
[reta, ceta] = size(DETA);
[retc, cetc] = size(DETC);
[rett, cett] = size(DETT);

% Total number of events in tag capture bands.
tcb = sum([cetg ceta cetc cett]);

% Display the tag capture statistics.
disp(['--> Number of all filtered events (exp): ', num2str(re)]);
disp(['--> Number of events in tag capture box (exp): ', num2str(u)]);
disp(['--> Number of events in tag capture bands (exp): ', num2str(tcb)]);
fprintf('\n');
disp(['--> Number of G-tag capture events in tag capture box (exp): ', num2str(cetg)]);
disp(['--> Number of A-tag capture events in tag capture box (exp): ', num2str(ceta)]);
disp(['--> Number of C-tag capture events in tag capture box (exp): ', num2str(cetc)]);
disp(['--> Number of T-tag capture events in tag capture box (exp): ', num2str(cett)]);
fprintf('\n');
disp(['--> Percent of G-tag capture events in tag capture box (exp): ', num2str(cetg/tcb*100)]);
disp(['--> Percent of A-tag capture events in tag capture box (exp): ', num2str(ceta/tcb*100)]);
disp(['--> Percent of C-tag capture events in tag capture box (exp): ', num2str(cetc/tcb*100)]);
disp(['--> Percent of T-tag capture events in tag capture box (exp): ', num2str(cett/tcb*100)]);

fprintf('\n');
disp('--> TAG CAPTURE ERROR ANALYSIS SECTION - BGR');

% Define array container for dwell time and current blockade level data of 
% 'background' data set = {G, A, C, T}.
DBTG = []; CBTG = [];
DBTA = []; CBTA = [];
DBTC = []; CBTC = [];
DBTT = []; CBTT = [];

% Counter of G/A/C/T-tag capture events.
gb = 1; ab = 1; bb = 1; tb = 1;

% Filter out events in the 'tag capture' for 'experiment' data set.
for n = 1:length(DB) 
    
    % If 'background event is in dwell time cutoff range (x) as well as in 
    % capture range (y), i.e., in tag capture box. 
    if (dwell_cutoff(1) < DB(n) && DB(n) < dwell_cutoff(2)) && ...
       (tag_capture(1) < CB(n) && CB(n) < tag_capture(8)) 

        % If event is in G-tag capture range.
        if tag_capture(1) < CB(n) && CB(n) < tag_capture(2) 
            
            % Remember these G-tag capture events.
            DBTG(gb) = DB(n);
            CBTG(gb) = CB(n);
            
            % Increment boxed-in G-tag capture event counter.
            gb = gb+1;
            
        end
        
        % If event is in A-tag capture range.
        if tag_capture(3) < CB(n) && CB(n) < tag_capture(4) 
            
            % Remember these A-tag capture events.
            DBTA(ab) = DB(n);
            CBTA(ab) = CB(n);
            
            % Increment boxed-in A-tag capture event counter.
            ab = ab+1;
            
        end
        
        % If event is in C-tag capture range.
        if tag_capture(5) < CB(n) && CB(n) < tag_capture(6) 
            
            % Remember these C-tag capture events.
            DBTC(bb) = DB(n);
            CBTC(bb) = CB(n);
            
            % Increment boxed-in C-tag capture event counter.
            bb = bb+1;
            
        end
        
        % If event is in T-tag capture range.
        if tag_capture(7) < CB(n) && CB(n) < tag_capture(8) 
            
            % Remember these T-tag capture events.
            DBTT(tb) = DB(n);
            CBTT(tb) = CB(n);
            
            % Increment boxed-in T-tag capture event counter.
            tb = tb+1;
            
        end
    end  
end

% Obtain size of all filtered and tag capture events.
[rb, cq] = size(DB);
[rbtg, cbtg] = size(DBTG);
[rbta, cbta] = size(DBTA);
[rbtc, cbtc] = size(DBTC);
[rbtt, cbtt] = size(DBTT);

% Display the tag capture statistics.
disp(['--> Number of all filtered events (bgr): ', num2str(rb)]);

disp(['--> Number of G-tag capture events in tag capture box (bgr): ', num2str(cbtg)]);
disp(['--> Number of A-tag capture events in tag capture box (bgr): ', num2str(cbta)]);
disp(['--> Number of C-tag capture events in tag capture box (bgr): ', num2str(cbtc)]);
disp(['--> Number of T-tag capture events in tag capture box (bgr): ', num2str(cbtt)]);
fprintf('\n');
disp(['--> Percent of G-tag capture events in tag capture box (bgr): ', num2str(cbtg/rb*100)]);
disp(['--> Percent of A-tag capture events in tag capture box (bgr): ', num2str(cbta/rb*100)]);
disp(['--> Percent of C-tag capture events in tag capture box (bgr): ', num2str(cbtc/rb*100)]);
disp(['--> Percent of T-tag capture events in tag capture box (bgr): ', num2str(cbtt/rb*100)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       EVENT ANALYSIS FINALIZATION                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Navigate to working directory.
cd(work_dir);

% Generate MAT-file structure for containing all these statistics.
m = matfile('err_stats_6-a4.mat', 'Writable', true);

% All experiment/background events.
m.dwellt_all_exp = DE;
m.dwellt_median_all_exp = median(DE);
m.dwellt_mean_all_exp = mean(DE);
m.dwellt_std_all_exp = std(DE);

m.ncurr_all_exp = CE;
m.ncurr_median_all_exp = median(CE);
m.ncurr_mean_all_exp = mean(CE);
m.ncurr_std_all_exp = std(CE);

m.dwellt_all_bgr = DB;
m.dwellt_median_all_bgr = median(DB);
m.dwellt_mean_all_bgr = mean(DB);
m.dwellt_std_all_bgr = std(DB);

m.ncurr_all_bgr = CB;
m.ncurr_median_all_bgr = median(CB);
m.ncurr_mean_all_bgr = mean(CB);
m.ncurr_std_all_bgr = std(CB);

% Tag capture zone experiment/background events only.
m.dwellt_tcz_exp_G = DETG';
m.ncurr_tcz_exp_G = CETG';
m.dwellt_tcz_bgr_G = DBTG';
m.ncurr_tcz_bgr_G = CBTG';

m.dwellt_tcz_exp_A = DETA';
m.ncurr_tcz_exp_A = CETA';
m.dwellt_tcz_bgr_A = DBTA';
m.ncurr_tcz_bgr_A = CBTA';

m.dwellt_tcz_exp_C = DETC';
m.ncurr_tcz_exp_C = CETC';
m.dwellt_tcz_bgr_C = DBTC';
m.ncurr_tcz_bgr_C = CBTC';

m.dwellt_tcz_exp_T = DETT';
m.ncurr_tcz_exp_T = CETT';
m.dwellt_tcz_bgr_T = DBTT';
m.ncurr_tcz_bgr_T = CBTT';

% Return the confusion matrix row.
E = ([ge ae cc te] - 1) / tcb * 100;
B = ([gb ab bb tb] - 1) / rb * 100;
A = [E; B];

fprintf('\n');
disp('--> Frequency analysis end');
fprintf('\n');
