%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: December 8, 2014.
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
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function hist_comp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           ANALYSIS STARTUP                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n');
disp('--> Histogram comparator start');
fprintf('\n');

% Set default number formatting.
format short;

% Define current working directory.
work_dir = pwd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     'EXPERIMENT' STATISTICS SECTION                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> EXPERIMENT STATISTICS SECTION');

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

for i = 1:length(list)
  
  disp(['--> Processing file: ', list(i).name]);  
  
  % Access data stack for dataset 'i'.
  stack = var.stack(i);	
  
  % Access data structure containing filter level 2 (all clean single 
  % channels) event information.
  cap = stack.stat(2).cap;

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

% Obtain size of all filtered events.
[r, c] = size(DE);

% Display the number of maximum events captured by a single pore.
disp(['--> Number of all filtered events: ', num2str(r)]);

% Create normalized current blockade level (pA) vs. dwell time (s) scatter 
% plot of all runs.
figure(1);
scatter(DE, CE, 'MarkerEdgeColor', 'red'); 
set(gca, 'XScale', 'log');
title('Experiment - Normalized Current Blockade Level vs. Dwell Time');
xlabel('Dwell Time (s)');
ylabel('Normalized Current Blockade Level (%)');
savefig('oc_vs_dwell_E.fig');
print('-dbmp', 'oc_vs_dwell_E.bmp'); 
disp('--> Scatter plot created for experiment dataset');
fprintf('\n');

% Create 2D plot of dwell time distribution with Gaussian fit. 
figure(2);
XE = logspace(-3, 3, 100);
hist(DE, XE);   
set(gca, 'XScale', 'log');
title('Experiment - Dwell Time Distribution');
xlabel('Dwell Time (s)');
ylabel('Count (#)');
xlim([1e-3 1e+3]);
h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black')
f = ezfit('gauss');
showfit(f, 'fitcolor', 'red', 'fitlinewidth', 2);
grid;
savefig('dwell_dist_E.fig');
print('-dbmp', 'dwell_dist_E.bmp'); 
fprintf('\n');
disp('--> Dwell time distribution 2D plot created for experiment dataset');
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                      'CONTROL' STATISTICS SECTION                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> CONTROL STATISTICS SECTION');

% Navigate to working directory.
cd(work_dir);

% Navigate to 'control' data directory.
if ~exist('control', 'dir')
  mkdir('control');
end

cd 'control';

% Load the event data structure into workspace.
load('var.mat');

% Read in all 'control' statistic text files one-by-one.
list = dir('data*_vardat');

% Define array container for dwell time data.
DC = []; 

% Define array container for current blockade level data.
CC = [];

for i = 1:length(list)
  
  disp(['--> Processing file: ', list(i).name]);  
  
  % Access data stack for dataset 'i'.
  stack = var.stack(i);	
  
  % Access data structure containing filter level 2 (all clean single 
  % channels) event information.
  cap = stack.stat(2).cap;

  % Retrieve all current blockade level and corresponding dwell time pairs.
  dwell = cap.dwell(:, :);
  imed = cap.imed(:, :);

  % Filter out NaNs and diplay 2D scatter plot and histogram with a 
  % distribution fit.
  iD = ~isnan(dwell);
  dc = dwell(iD(:));
  cc = imed(iD(:));
  
  % Obtain size of maximum events.
  [r, c] = size(dc);
  
  % Display the number of maximum events captured by a single pore.
  disp(['--> Maximum number of events: ', num2str(r)]);
  fprintf('\n');
    
  % Collect all filtered events and update storage array.
  DC = [DC; dc];
  CC = [CC; cc];
  
end

% Obtain size of all filtered events.
[r, c] = size(DC);

% Display the number of maximum events captured by a single pore.
disp(['--> Number of all filtered events: ', num2str(r)]);

% Create normalized current blockade level (pA) vs. dwell time (s) scatter 
% plot of all runs.
figure(3);
scatter(DC, CC, 'MarkerEdgeColor', 'blue'); 
set(gca, 'XScale', 'log');
title('Control - Normalized Current Blockade Level vs. Dwell Time');
xlabel('Dwell Time (s)');
ylabel('Normalized Current Blockade Level (%)');
savefig('oc_vs_dwell_C.fig');
print('-dbmp', 'oc_vs_dwell_C.bmp'); 
disp('--> Scatter plot created for control dataset');
fprintf('\n');

% Create 2D plot of dwell time distribution with Gaussian fit. 
figure(4);
XC = logspace(-3, 3, 100);
hist(DC, XC);   
set(gca, 'XScale', 'log');
title('Control - Dwell Time Distribution');
xlabel('Dwell Time (s)');
ylabel('Count (#)');
xlim([1e-3 1e+3]);
h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [0.7 0.7 1], 'EdgeColor', 'black');
f = ezfit('gauss');
showfit(f, 'fitcolor', 'blue', 'fitlinewidth', 2);
grid;
savefig('dwell_dist_C.fig');
print('-dbmp', 'dwell_dist_C.bmp'); 
fprintf('\n');
disp('--> Dwell time distribution 2D plot created for control dataset');
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                      'OVERLAY' STATISTICS SECTION                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> OVERLAY STATISTICS SECTION');

% Navigate to working directory.
cd(work_dir);

% Create normalized current blockade level (pA) vs. dwell time (s) scatter 
% plot of all runs overlaying 'control' with 'experiment' data.
figure(5);
scatter(DC, CC, 'MarkerEdgeColor', 'blue'); 
hold on;
scatter(DE, CE, 'MarkerEdgeColor', 'red'); 
set(gca, 'XScale', 'log');
title('Overlay - Normalized Current Blockade Level vs. Dwell Time');
xlabel('Dwell Time (s)');
ylabel('Normalized Current Blockade Level (%)');
savefig('oc_vs_dwell_O.fig');
print('-dbmp', 'oc_vs_dwell_O.bmp'); 
disp('--> Scatter plot created for overlaid dataset');

% Create 2D plot of dwell time distribution with Gaussian fit overlaying 
% 'control' with 'experiment' data.
figure(6);
hist(DC, XC);   
hold on;
hist(DE, XE);
h = findobj(gca, 'Type', 'patch');
set(h(1), 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black');
set(h(2), 'FaceColor', [0.7 0.7 1], 'EdgeColor', 'black');
set(gca, 'XScale', 'log');
grid;
title('Overlay - Dwell Time Distribution');
xlabel('Dwell Time (s)');
ylabel('Count (#)');
xlim([1e-3 1e+3]);
savefig('dwell_dist_O.fig');
print('-dbmp', 'dwell_dist_O.bmp'); 
disp('--> Dwell time distribution 2D plot created for overlaid dataset');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          ANALYSIS FINALIZATION                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n');
disp('--> Histogram comparator end');
fprintf('\n');
