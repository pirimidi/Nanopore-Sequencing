%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: August 24, 2015.
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: Given 4 'stats.mat' structure of static capture vs. background
% data sets, this program iterates through all of them and generates 8 box-
% plots for:
%
% (1) DE := dwell time of all filtered events (exp)
% (2) CE := normalized current of all filtered events (exp)
% (3) DB := dwell time of all filtered events (bgr)
% (4) CB := normalized current of all filtered events (bgr)
% (5) DET := dwell time of events in tag cature zone (exp)
% (6) CET := normalized current of events in tag cature zone (exp)
% (7) DBT := dwell time of events in tag capture zone (bgr)
% (8) CBT := normalized current of events in tag capture zone (bgr) 
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function boxplot_svb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         EVENT ANALYSIS STARTUP                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all')

fprintf('\n');
disp('--> Boxplots for static capture versus background start');
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
%                       BOXPLOT STATISTICS SECTION                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> BOXPLOT STATISTICS SECTION');

% Navigate to 'stats' data directory.
if ~exist('stats', 'dir')
  mkdir('stats');
end

cd 'stats';

% Read in all 'stats' text file names one-by-one.
list = dir('stats_*');

% Define array container for dwell time, normalized current and size of all 
% filtered events (bgr).
DB = []; 
CB = [];
SB = [];

% Define array container for dwell time, normalized current and size of 
% events in tag cature zone (bgr).
DBT = [];
CBT = [];
SBT = [];

for i = 1:length(list)

  disp(['--> Processing structure: ', list(i).name]);  

  % Load the stats data structure into workspace.
  load(list(i).name);

  % Generate container arrays for all filtered events (bgr).
  DB = [DB; dwellt_all_bgr];
  CB = [CB; ncurr_all_bgr];
  
  % Determine number of events in this data set.
  SB(i) = length(dwellt_all_bgr);    
  
  % Generate container arrays for events in tag capture zone (bgr).
  DBT = [DBT; dwellt_tcz_bgr];
  CBT = [CBT; ncurr_tcz_bgr];  
  
  % Determine number of events in this data set.
  SBT(i) = length(dwellt_tcz_bgr);    

end

% Build grouping variable array.
bases = ['G', 'A', 'C', 'T'];

GB = [];
GBT = [];

for j = 1:length(list) 
    
    % Initialize containers.
    gb = [];
    gbt = [];
    
    for k = 1:SB(j)
        gb{k} = bases(j);
    end
    
    for k = 1:SBT(j)
        gbt{k} = bases(j);
    end
    
    % Update concatinated grouping variable array.
    GB = [GB; gb'];
    GBT = [GBT; gbt'];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        BOXPLOT ANALYSIS FINALIZATION                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> GENERATING BOXPLOTS');

% Navigate to plotting directory.
cd '../plots';

% Create dwell time (s) boxplot of all filtered events (bgr).
figure(1);
boxplot(DB, GB); 
set(gca, 'YScale', 'log', 'ylim', [1E-3 1E+3]);
title('Background - Dwell Time - All Filtered Events');
xlabel('Tagged Nucleotides');
ylabel('Dwell Time (s)');
savefig('box_dwell_all_bgr.fig');
print('-dbmp', 'box_dwell_all_bgr.bmp'); 
disp('--> Boxplot created for all filtered events (bgr)');

% Create normalized current (%) boxplot of all filtered events (bgr).
figure(2);
boxplot(CB, GB); 
set(gca, 'ylim', [0.0 0.7001]);
title('Background - Normalized Current Blockade Level (%) - All Filtered Events');
xlabel('Tagged Nucleotides');
ylabel('Normalized Current Blockade Level (%)');
savefig('box_ncurr_all_bgr.fig');
print('-dbmp', 'box_ncurr_all_bgr.bmp'); 
disp('--> Boxplot created for all filtered events (bgr)');

%%%

% Create dwell time (s) boxplot of tag capture zone events (bgr).
figure(3);
boxplot(DBT, GBT); 
set(gca, 'YScale', 'log', 'ylim', [1E-3 1E+3]);
title('Background - Dwell Time - Tag Capture Zone');
xlabel('Tagged Nucleotides');
ylabel('Dwell Time (s)');
savefig('box_dwell_tcz_bgr.fig');
print('-dbmp', 'box_dwell_tcz_bgr.bmp'); 
disp('--> Boxplot created for tag capture zone events (bgr)');

% Create normalized current (%) boxplot of tag capture zone events (bgr).
figure(4);
boxplot(CBT, GBT);
set(gca, 'ylim', [0.0 0.7001]);
title('Background - Normalized Current Blockade Level (%) - Tag Capture Zone');
xlabel('Tagged Nucleotides');
ylabel('Normalized Current Blockade Level (%)');
savefig('box_ncurr_tcz_bgr.fig');
print('-dbmp', 'box_ncurr_tcz_bgr.bmp'); 
disp('--> Boxplot created for tag capture zone events (bgr)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        BOXPLOT ANALYSIS FINALIZATION                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> Boxplots for static capture versus background end');
fprintf('\n');
