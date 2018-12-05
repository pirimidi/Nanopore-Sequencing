%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: August 27, 2015.
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
% "How many times are we in any of the 3 wrong signal band(s) for all 
% capture events?" (Answer is given as a % of total events.)
%
% Input arguments:
%
% (1) filter_e := filter level for pores in experimrntal data [1-4]
% (2) order_e := stack order, i.e., experimental data sets to access
% (3) tag_capture := lower and upper tag capture limits, an array of {GACT}
% (4) dwell_cutoff := background dwell time cutoff (in sec), an array
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function [E] = poly_error(filter_e, order_e, tag_capture, dwell_cutoff)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         ERROR ANALYSIS STARTUP                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all')

fprintf('\n');
disp('--> Polymerase incoporation error analysis start');
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

% Obtain size of all filtered events.
[r, c] = size(DE);

% Display the number of maximum events captured by a single pore.
disp(['--> Number of all filtered events (exp): ', num2str(r)]);
fprintf('\n');

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
for m = 1:length(DE)

    % If event is in dwell time cutoff range (x) as well as in capture 
    % range (y), i.e., in tag capture box. 
    if (dwell_cutoff(1) < DE(m) && DE(m) < dwell_cutoff(2)) && ...
       (tag_capture(1) < CE(m) && CE(m) < tag_capture(8))  
   
        % Increment events in tag capture box.
        u = u+1;
    
        % If event is in G-tag capture range.
        if tag_capture(1) < CE(m) && CE(m) < tag_capture(2) 
            
            % Remember these G-tag capture events.
            DETG(ge) = DE(m);
            CETG(ge) = CE(m);
            
            % Increment boxed-in G-tag capture event counter.
            ge = ge+1;
            
        end
        
        % If event is in A-tag capture range.
        if tag_capture(3) < CE(m) && CE(m) < tag_capture(4) 
            
            % Remember these A-tag capture events.
            DETA(ae) = DE(m);
            CETA(ae) = CE(m);
            
            % Increment boxed-in A-tag capture event counter.
            ae = ae+1;
            
        end
        
        % If event is in C-tag capture range.
        if tag_capture(5) < CE(m) && CE(m) < tag_capture(6) 
            
            % Remember these C-tag capture events.
            DETC(cc) = DE(m);
            CETC(cc) = CE(m);
            
            % Increment boxed-in C-tag capture event counter.
            cc = cc+1;
            
        end
        
        % If event is in T-tag capture range.
        if tag_capture(7) < CE(m) && CE(m) < tag_capture(8) 
            
            % Remember these T-tag capture events.
            DETT(te) = DE(m);
            CETT(te) = CE(m);
            
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       EVENT ANALYSIS FINALIZATION                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Navigate to working directory.
cd(work_dir);

% Generate MAT-file structure for containing all these statistics.
m = matfile('err_stats_6-d3.mat', 'Writable', true);

% All experiment/background events.
m.dwellt_all_exp = DE;
m.dwellt_median_all_exp = median(DE);
m.dwellt_mean_all_exp = mean(DE);
m.dwellt_std_all_exp = std(DE);

m.ncurr_all_exp = CE;
m.ncurr_median_all_exp = median(CE);
m.ncurr_mean_all_exp = mean(CE);
m.ncurr_std_all_exp = std(CE);

% Tag capture zone experiment/background events only.
m.dwellt_tcz_exp_G = DETG';
m.ncurr_tcz_exp_G = CETG';

m.dwellt_tcz_exp_A = DETA';
m.ncurr_tcz_exp_A = CETA';

m.dwellt_tcz_exp_C = DETC';
m.ncurr_tcz_exp_C = CETC';

m.dwellt_tcz_exp_T = DETT';
m.ncurr_tcz_exp_T = CETT';

% Return the confusion matrix row.
E = ([ge ae cc te] - 1) / tcb * 100;

fprintf('\n');
disp('--> Polymerase incoporation error analysis end');
fprintf('\n');
