%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: August 18, 2015.
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: Given a 'var.mat' structure of a 'vardat' run for N data sets,
% this program iterates through all of them and extracts filter level 2 
% (all clean pores) event information of current blockade level and 
% corresponding dwell time pairs, then generates a 2D scatter plot and 
% histogram with a distribution fit.
%
% Input arguments:
%
% (1) growth_b := scaling factor for background data
% (2) filter_b := filter level for pores in background data [1-4]
% (3) order_b := stack order, i.e., background data sets to access
% (4) tag_capture := lower and upper tag capture limits, an array
% (5) dwell_cutoff := background dwell time cutoff (in sec) 
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function tag_bgr(growth_b, filter_b, order_b, tag_capture, dwell_cutoff)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         EVENT ANALYSIS STARTUP                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all')

fprintf('\n');
disp('--> Tag background start');
fprintf('\n');

% Set default number formatting.
format short;

% Define current working directory.
work_dir = pwd;

% Navigate to 'plots' directory.
if ~exist('plots', 'dir')
  mkdir('plots');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  'BACKGROUND' EVENTS STATISTICS SECTION                 %
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

for k = 1:growth_b
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
end

% Obtain size of all filtered events.
[r, c] = size(DB);

% Display the number of maximum events captured by a single pore.
disp(['--> Number of all filtered events: ', num2str(r)]);

% Navigate to plotting directory.
cd '../plots';

% Create normalized current blockade level (pA) vs. dwell time (s) scatter 
% plot of all runs.
figure(1);
scatter(DB, CB, 'MarkerEdgeColor', 'blue'); 
set(gca, 'XScale', 'log');
title('Background - Normalized Current Blockade Level vs. Dwell Time');
xlabel('Dwell Time (s)');
ylabel('Normalized Current Blockade Level (%)');
savefig('oc_vs_dwell_B.fig');
print('-dbmp', 'oc_vs_dwell_B.bmp'); 
disp('--> Scatter plot created for background dataset');

% Create 2D plot of dwell time distribution with Gaussian fit. 
figure(2);
XB = logspace(-3, 3, 100);
hist(DB, XB);   
set(gca, 'XScale', 'log');
title('Background - Dwell Time Distribution');
xlabel('Dwell Time (s)');
ylabel('Count (#)');
xlim([1e-3 1e+3]);
h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [0.7 0.7 1], 'EdgeColor', 'black');
% f = ezfit('gauss');
% showfit(f, 'fitcolor', 'blue', 'fitlinewidth', 2);
grid;
savefig('dwell_dist_B.fig');
print('-dbmp', 'dwell_dist_B.bmp'); 
disp('--> Dwell time distribution 2D plot created for background dataset');
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  'TAG CAPTURE' EVENT STATISTICS SECTION                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> TAG CAPTURE EVENT STATISTICS SECTION');

% Define array container for dwell time and current blockade level data of 
% 'background' data set.
DBT = []; 
CBT = [];

% Counter of events > dwell time cutoff.
v = 0;

% Counter of tag raw_template events.
g = 1;

% Filter out events in the 'tag capture' for 'background' data set.
for n = 1:length(DB)
    
    % If event is less than dwell time cutoff. 
    if DB(n) > dwell_cutoff

        % Increment events > dwell time cutoff counter.
        v = v+1;        
    end    

    % If event is in tag capture.
    if tag_capture(1) < CB(n) && CB(n) < tag_capture(2) 
        
        % If event is less than dwell time cutoff. 
        if DB(n) > dwell_cutoff
            
            % Remember these tag capture events.
            DBT(g) = DB(n);
            CBT(g) = CB(n);
            
            % Increment tag capture event counter.
            g = g+1;
            
        end
    end  
end

% Obtain size of all filtered and tag capture events.
[ra, ca] = size(DB);
[rt, ct] = size(DBT);

% Display the tag capture statistics.
disp(['--> Number of all filtered events (bgr): ', num2str(ra)]);
disp(['--> Number of events > dwell time cutoff (bgr): ', num2str(v)]);
disp(['--> Number of tag capture events > dwell time cutoff (bgr): ', num2str(ct)]);
disp(['--> Percent of tag capture events to all filtered events (bgr): ', num2str(ct/ra*100)]);
disp(['--> Percent of tag capture events to events > dwell time cutoff (bgr): ', num2str(ct/v*100)]);
fprintf('\n');

% Create normalized current blockade level (pA) vs. dwell time (s) scatter 
% plot of all runs overlaying 'background' with 'experiment' data.
figure(3);
scatter(DBT, CBT, 'MarkerEdgeColor', 'blue');
set(gca, 'XScale', 'log');
title('Tag capture - Normalized Current Blockade Level vs. Dwell Time');
xlabel('Dwell Time (s)');
ylabel('Normalized Current Blockade Level (%)');
savefig('oc_vs_dwell_T.fig');
print('-dbmp', 'oc_vs_dwell_T.bmp'); 
disp('--> Scatter plot created for tag capture dataset');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       EVENT ANALYSIS FINALIZATION                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> Static capture to background event comparator end');
fprintf('\n');
