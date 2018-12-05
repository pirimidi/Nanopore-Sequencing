%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: August 20, 2015.
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
% (1) growth_e := scaling factor for experimental data
% (2) growth_c := scaling factor for background data
% (3) filter_e := filter level for pores in experimrntal data [1-4]
% (4) filter_c := filter level for pores in background data [1-4]
% (5) order_e := stack order, i.e., experimental data sets to access
% (6) order_c := stack order, i.e., background data sets to access
% (7) tag_capture := lower and upper tag capture limits, an array
% (8) dwell_cutoff := background dwell time cutoff (in sec) 
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function statcap_vs_bgr(growth_e, growth_b, filter_e, filter_b, ...
                        order_e, order_b, tag_capture, dwell_cutoff)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         EVENT ANALYSIS STARTUP                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all')

fprintf('\n');
disp('--> Static capture to background event comparator start');
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
%                  'EXPERIMENT' EVENT STATISTICS SECTION                  %
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

for i = 1:growth_e
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
end

% Obtain size of all filtered events.
[r, c] = size(DE);

% Display the number of maximum events captured by a single pore.
disp(['--> Number of all filtered events: ', num2str(r)]);

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

% Navigate to plotting directory.
cd '../plots';

% Create normalized current blockade level (pA) vs. dwell time (s) scatter 
% plot of all runs.
figure(1);
scatter(DEU, CEU, 'MarkerEdgeColor', 'red'); 
set(gca, 'XScale', 'log');
axis([1e-3 1e+3 0 0.7001]);
title('Experiment - Normalized Current Blockade Level vs. Dwell Time');
xlabel('Dwell Time (s)');
ylabel('Normalized Current Blockade Level (%)');
savefig('oc_vs_dwell_E.fig');
print('-dbmp', 'oc_vs_dwell_E.bmp'); 
disp('--> Scatter plot created for experiment dataset');

% Create histogram of dwell time with Gaussian fit. 
figure(2);
XE = logspace(-3, 3, 100);
hist(DEU, XE);   
set(gca, 'XScale', 'log');
title('Experiment - Dwell Time Distribution');
xlabel('Dwell Time (s)');
ylabel('Count (#)');
xlim([1e-3 1e+3]);
h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black')
% f = ezfit('gauss');
% showfit(f, 'fitcolor', 'red', 'fitlinewidth', 2);
grid;
savefig('dwell_dist_E.fig');
print('-dbmp', 'dwell_dist_E.bmp'); 
disp('--> Dwell time histogram created for experiment dataset');

% Create histogram of normalized current blockade level with Gaussian fit. 
figure(3);
hist(CEU, 100);   
title('Experiment - Normalized Current Blockade Level Distribution');
xlabel('Normalized Current Blockade Level (%)');
ylabel('Count (#)');
%xlim([0 0.7]);
h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black')
% f = ezfit('gauss');
% showfit(f, 'fitcolor', 'red', 'fitlinewidth', 2);
grid;
savefig('ncurr_dist_E.fig');
print('-dbmp', 'ncurr_dist_E.bmp'); 
disp('--> Normalized current historgram created for experiment dataset');
fprintf('\n');

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
figure(4);
scatter(DB, CB, 'MarkerEdgeColor', 'blue'); 
set(gca, 'XScale', 'log');
axis([1e-3 1e+3 0 0.7001]);
title('Background - Normalized Current Blockade Level vs. Dwell Time');
xlabel('Dwell Time (s)');
ylabel('Normalized Current Blockade Level (%)');
savefig('oc_vs_dwell_B.fig');
print('-dbmp', 'oc_vs_dwell_B.bmp'); 
disp('--> Scatter plot created for background dataset');

% Create histogram of dwell time with Gaussian fit. 
figure(5);
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
disp('--> Dwell time histogram created for background dataset');

% Create histogram of normalized current blockade level with Gaussian fit. 
figure(6);
hist(CB, 100);   
title('Background - Normalized Current Blockade Level Distribution');
xlabel('Normalized Current Blockade Level (%)');
ylabel('Count (#)');
xlim([0 0.7]);
h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [0.7 0.7 1], 'EdgeColor', 'black');
f = ezfit('gauss');
showfit(f, 'fitcolor', 'blue', 'fitlinewidth', 2);
grid;
savefig('ncurr_dist_B.fig');
print('-dbmp', 'ncurr_dist_B.bmp'); 
fprintf('\n');
disp('--> Normalized current histogram created for background dataset');
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    'OVERLAY' EVENT STATISTICS SECTION                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> OVERLAY STATISTICS SECTION');

% Navigate to plotting directory.
cd '../plots';

% Create normalized current blockade level (pA) vs. dwell time (s) scatter 
% plot of all runs overlaying 'background' with 'experiment' data.
figure(7);
scatter(DB, CB, 'MarkerEdgeColor', 'blue');
hold on;
scatter(DEU, CEU, 'MarkerEdgeColor', 'red');
set(gca, 'XScale', 'log');
axis([1e-3 1e+3 0 0.7001]);
title('Overlay - Normalized Current Blockade Level vs. Dwell Time');
xlabel('Dwell Time (s)');
ylabel('Normalized Current Blockade Level (%)');
savefig('oc_vs_dwell_O.fig');
print('-dbmp', 'oc_vs_dwell_O.bmp'); 
disp('--> Scatter plot created for overlaid dataset');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     'TAG CAPTURE' EVENT STATISTICS SECTION              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> TAG CAPTURE EVENT STATISTICS SECTION');

% Define array container for dwell time and current blockade level data of 
% 'experiment' data set.
DET = []; 
CET = [];

% Counter of events > dwell time cutoff.
u = 0;

% Counter of tag capture events.
h = 1;

% Filter out events in the 'tag capture' for 'experiment' data set.
for m = 1:length(DE)
    
    % If event is less than dwell time cutoff. 
    if DE(m) > dwell_cutoff

        % Increment events > dwell time cutoff counter.
        u = u+1;        
    end

    % If event is in tag capture.
    if tag_capture(1) < CE(m) && CE(m) < tag_capture(2) 
        
        % If event is less than dwell time cutoff. 
        if DE(m) > dwell_cutoff
            
            % Remember these tag capture events.
            DET(h) = DE(m);
            CET(h) = CE(m);
            
            % Increment tag capture event counter.
            h = h+1;
            
        end
    end  
end

% Obtain size of all filtered and tag capture events.
[re, ce] = size(DE);
[ret, cet] = size(DET);

% Display the tag capture statistics.
disp(['--> Number of all filtered events (exp): ', num2str(re)]);
disp(['--> Number of events > dwell time cutoff (exp): ', num2str(u)]);
disp(['--> Number of tag capture events > dwell time cutoff (exp): ', num2str(cet)]);
disp(['--> Percent of tag capture events to events > dwell time cutoff (exp): ', num2str(cet/u*100)]);
disp(['--> Percent of tag capture events to total (exp): ', num2str(cet/re*100)]);
fprintf('\n');

% Define array container for dwell time and current blockade level data of 
% 'experiment' data set.
DBT = []; 
CBT = [];

% Counter of events > dwell time cutoff.
v = 0;

% Counter of tag raw_template events.
g = 1;

% Filter out events in the 'tag capture' for 'experiment' data set.
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
[rb, cb] = size(DB);
[rbt, cbt] = size(DBT);

% Display the tag capture statistics.
disp(['--> Number of all filtered events (bgr): ', num2str(rb)]);
disp(['--> Number of events > dwell time cutoff (bgr): ', num2str(v)]);
disp(['--> Number of tag capture events > dwell time cutoff (bgr): ', num2str(cbt)]);
disp(['--> Percent of tag capture events to events > dwell time cutoff (bgr): ', num2str(cbt/v*100)]);
disp(['--> Percent of tag capture events to all filtered events (bgr): ', num2str(cbt/rb*100)]);
fprintf('\n');

% Navigate to plotting directory.
cd '../plots';

% Create normalized current blockade level (pA) vs. dwell time (s) scatter 
% plot of all runs overlaying 'background' with 'experiment' data.
figure(8);
scatter(DET, CET, 'MarkerEdgeColor', 'red'); 
hold on;
scatter(DBT, CBT, 'MarkerEdgeColor', 'blue');
set(gca, 'XScale', 'log');
%axis([1e-3 1e+3 0 0.7001]);
title('Tag capture - Normalized Current Blockade Level vs. Dwell Time');
xlabel('Dwell Time (s)');
ylabel('Normalized Current Blockade Level (%)');
savefig('oc_vs_dwell_T.fig');
print('-dbmp', 'oc_vs_dwell_T.bmp'); 
disp('--> Scatter plot created for tag capture dataset');

% Create histogram of dwell time with Gaussian fit overlaying 'background' 
% with 'experiment' data.
figure(9);
hist(DET, XE);
hold on;
hist(DBT, XB);
h = findobj(gca, 'Type', 'patch');
set(h(1), 'FaceColor', [0.7 0.7 1], 'EdgeColor', 'black');
set(h(2), 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black');
set(gca, 'XScale', 'log');
grid;
title('Tag capture - Dwell Time Distribution');
xlabel('Dwell Time (s)');
ylabel('Count (#)');
xlim([1e-3 1e+3]);
savefig('dwell_dist_T.fig');
print('-dbmp', 'dwell_dist_T.bmp'); 
disp('--> Dwell time histogram created for tag capture dataset');

% Create histogram of normalized current blockade level with Gaussian fit 
% overlaying 'background' with 'experiment' data.
figure(10);
hist(CET, 100);
hold on;
hist(CBT, 100);
h = findobj(gca, 'Type', 'patch');
set(h(1), 'FaceColor', [0.7 0.7 1], 'EdgeColor', 'black');
set(h(2), 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black');
grid;
title('Tag capture - Normalized Current Blockade Level Distribution');
xlabel('Normalized Current Blockade Level (%)');
ylabel('Count (#)');
%xlim([0 0.7]);
savefig('ncurr_dist_T.fig');
print('-dbmp', 'ncurr_dist_T.bmp'); 
disp('--> Normalized current histogram for tag capture dataset');
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       EVENT ANALYSIS FINALIZATION                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Navigate to working directory.
cd(work_dir);

disp('--> EVENT ANALYSIS FINALIZATION');

% Display the median, mean and standard deviation of dwell time and 
% Normalized current for 'experiment' events.
disp(['--> Dwell time median of all filtered events (exp): ', num2str(median(DE))]);
disp(['--> Dwell time mean of all filtered events (exp): ', num2str(mean(DE))]);
disp(['--> Dwell time standard deviation of all filtered events (exp): ', num2str(std(DE))]);
fprintf('\n');
disp(['--> Normalized current median of all filtered events (exp): ', num2str(median(CE))]);
disp(['--> Normalized current mean of all filtered events (exp): ', num2str(mean(CE))]);
disp(['--> Normalized current standard deviation of all filtered events (exp): ', num2str(std(CE))]);
fprintf('\n');

% Display the median, mean and standard deviation of dwell time and 
% Normalized current for 'background' events.
disp(['--> Dwell time median of all filtered events (bgr): ', num2str(median(DB))]);
disp(['--> Dwell time mean of all filtered events (bgr): ', num2str(mean(DB))]);
disp(['--> Dwell time standard deviation of all filtered events (bgr): ', num2str(std(DB))]);
fprintf('\n');
disp(['--> Normalized current median of all filtered events (bgr): ', num2str(median(CB))]);
disp(['--> Normalized current mean of all filtered events (bgr): ', num2str(mean(CB))]);
disp(['--> Normalized current standard deviation of all filtered events (bgr): ', num2str(std(CB))]);
fprintf('\n');

% Display the median, mean and standard deviation of dwell time and 
% Normalized current for 'experiment' events in tag capture zone.
disp(['--> Dwell time median of events in tag cature zone (exp): ', num2str(median(DET))]);
disp(['--> Dwell time mean of events in tag cature zone (exp): ', num2str(mean(DET))]);
disp(['--> Dwell time standard deviation of events in tag capture zone (exp): ', num2str(std(DET))]);
fprintf('\n');
disp(['--> Normalized current median of events in tag cature zone (exp): ', num2str(median(CET))]);
disp(['--> Normalized current mean of events in tag cature zone (exp): ', num2str(mean(CET))]);
disp(['--> Normalized current standard deviation of events in tag capture zone (exp): ', num2str(std(CET))]);
fprintf('\n');

% Display the median, mean and standard deviation of dwell time and 
% Normalized current for 'background' events in tag capture zone.
disp(['--> Dwell time median of events in tag capture zone (bgr): ', num2str(median(DBT))]);
disp(['--> Dwell time mean of events in tag capture zone (bgr): ', num2str(mean(DBT))]);
disp(['--> Dwell time standard deviation of events in tag capture zone (bgr): ', num2str(std(DBT))]);
fprintf('\n');
disp(['--> Normalized current median of events in tag capture zone (bgr): ', num2str(median(CBT))]);
disp(['--> Normalized current mean of events in tag capture zone (bgr): ', num2str(mean(CBT))]);
disp(['--> Normalized current standard deviation of events in tag capture zone (bgr): ', num2str(std(CBT))]);

% Generate MAT-file structure for containing all these statistics.
m = matfile('stats_6-4.mat', 'Writable', true);

m.all_events_exp = re;
m.dwellt_events_exp = u;
m.tagcap_events_exp = cet;
m.tagcap_to_dwellt_exp = cet/u*100;
m.tagcap_to_all_exp = cet/re*100;

m.all_events_bgr = rb;
m.dwellt_events_bgr = v;
m.tagcap_events_bgr = cbt;
m.tagcap_to_dwellt_bgr = cbt/v*100;
m.tagcap_to_all_bgr = cbt/rb*100;

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
m.dwellt_tcz_exp = DET';
m.dwellt_median_tcz_exp = median(DET);
m.dwellt_mean_tcz_exp = mean(DET);
m.dwellt_std_tcz_exp = std(DET);

m.ncurr_tcz_exp = CET';
m.ncurr_median_tcz_exp = median(CET);
m.ncurr_mean_tcz_exp = mean(CET);
m.ncurr_std_tcz_exp = std(CET);

m.dwellt_tcz_bgr = DBT';
m.dwellt_median_tcz_bgr = median(DBT);
m.dwellt_mean_tcz_bgr = mean(DBT);
m.dwellt_std_tcz_bgr = std(DBT);

m.ncurr_tcz_bgr = CBT';
m.ncurr_median_tcz_bgr = median(CBT);
m.ncurr_mean_tcz_bgr = mean(CBT);
m.ncurr_std_tcz_bgr = std(CBT);

fprintf('\n');
disp('--> Static capture to background event comparator end');
fprintf('\n');
