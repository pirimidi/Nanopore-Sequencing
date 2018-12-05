%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: September 11, 2015.
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: Given a 'stats.mat' structure of 2-tag transition data sets, 
% this program iterates through all of them and generates 3 + 16 box-plots 
% for:
%
% (1) all filtered events
% (2) N->N transitions 
% (3) N->N+1 transitions 
% (4)-(19) pairwise combination of all {G, A, C, T} transitions
%
% for three parameters:
%
% (a) dwell time
% (b) wait time 
% (c) normalized current 
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function boxplot_transition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        STUTTER BOXPLOT STARTUP                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all')

fprintf('\n');
disp('--> Boxplots for tag transition analysis start');
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
%                    STUTTER BOXPLOT STATISTICS SECTION                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> STUTTER BOXPLOT STATISTICS SECTION');

% Navigate to 'stats' data directory.
if ~exist('stats', 'dir')
  mkdir('stats');
end

cd 'stats';

% Read in all 'stats' text file names one-by-one.
data = dir('*.mat'); 
load(data.name);

% Define array container for dwell time, wait time and normalized current
% of transitions.
D = []; 
W = [];
S = [];

% Define array container for dwell time, wait time of transitions.
D = [dwellt_all; dwellt_NN'; dwellt_NM'; ...
     dwellt_GG'; dwellt_GA'; dwellt_GC'; dwellt_GT'; ...
     dwellt_AG'; dwellt_AA'; dwellt_AC'; dwellt_AT'; ...
     dwellt_CG'; dwellt_CA'; dwellt_CC'; dwellt_CT'; ...
     dwellt_TG'; dwellt_TA'; dwellt_TC'; dwellt_TT'];
 
W = [wait_all; wait_NN'; wait_NM'; ...
     wait_GG'; wait_GA'; wait_GC'; wait_GT'; ...
     wait_AG'; wait_AA'; wait_AC'; wait_AT'; ...
     wait_CG'; wait_CA'; wait_CC'; wait_CT'; ...
     wait_TG'; wait_TA'; wait_TC'; wait_TT'];
 
C = [ncurr_all; ncurr_NN'; ncurr_NM'; ...
     ncurr_GG'; ncurr_GA'; ncurr_GC'; ncurr_GT'; ...
     ncurr_AG'; ncurr_AA'; ncurr_AC'; ncurr_AT'; ...
     ncurr_CG'; ncurr_CA'; ncurr_CC'; ncurr_CT'; ...
     ncurr_TG'; ncurr_TA'; ncurr_TC'; ncurr_TT'];
 
% Define array container for number of dwell time transitions.
LD = [length(dwellt_all); length(dwellt_NN); length(dwellt_NM); ...
      length(dwellt_GG); length(dwellt_GA); length(dwellt_GC); length(dwellt_GT); ...
      length(dwellt_AG); length(dwellt_AA); length(dwellt_AC); length(dwellt_AT); ...
      length(dwellt_CG); length(dwellt_CA); length(dwellt_CC); length(dwellt_CT); ...
      length(dwellt_TG); length(dwellt_TA); length(dwellt_TC); length(dwellt_TT)];
  
LW = [length(wait_all); length(wait_NN); length(wait_NM); ...
      length(wait_GG); length(wait_GA); length(wait_GC); length(wait_GT); ...
      length(wait_AG); length(wait_AA); length(wait_AC); length(wait_AT); ...
      length(wait_CG); length(wait_CA); length(wait_CC); length(wait_CT); ...
      length(wait_TG); length(wait_TA); length(wait_TC); length(wait_TT)];
  
LC = [length(ncurr_all); length(ncurr_NN); length(ncurr_NM); ...
      length(ncurr_GG); length(ncurr_GA); length(ncurr_GC); length(ncurr_GT); ...
      length(ncurr_AG); length(ncurr_AA); length(ncurr_AC); length(ncurr_AT); ...
      length(ncurr_CG); length(ncurr_CA); length(ncurr_CC); length(ncurr_CT); ...
      length(ncurr_TG); length(ncurr_TA); length(ncurr_TC); length(ncurr_TT)];

% Build grouping variable array.
trans = {'all', 'N->N', 'N->M', ...
         'G->G', 'G->A', 'G->C', 'G->T', ...
         'A->G', 'A->A', 'A->C', 'A->T', ...
         'C->G', 'C->A', 'C->C', 'C->T', ...
         'T->G', 'T->A', 'T->C', 'T->T'};
     
% Initialize grouping variable arrays.
DG = []; WG = []; CG = [];

for j = 1:length(LD)

    % Initialize containers.
    dg = {}; wg = {}; cg = {};
    
    for k = 1:LD(j)
        dg{k} = trans{j};
    end
    
    for k = 1:LW(j)
        wg{k} = trans{j};
    end

    for k = 1:LC(j)
        cg{k} = trans{j};
    end

    % Update concatinated grouping variable array.
    DG = [DG; dg'];
    WG = [WG; wg'];
    CG = [CG; cg'];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  STUTTER BOXPLOT ANALYSIS FINALIZATION                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> GENERATING STUTTER BOXPLOTS');

% Navigate to plotting directory.
cd '../plots';

% Create dwell time (s) boxplot of all transitions.
figure(1);
boxplot(D, DG); 
set(gca, 'ylim', [-0.1 0.8]);
%set(gca, 'YScale', 'log', 'ylim', [1E-3 1E+3]);
title('Dwell Time Transitions');
xlabel('Transitions');
ylabel('Dwell Time (s)');
savefig('box_dwell.fig');
print('-dbmp', 'box_dwell.bmp'); 
disp('--> Boxplot created for dwell time transitions');

% Create wait time (s) boxplot of all transitions.
figure(2);
boxplot(W, WG); 
set(gca, 'ylim', [-0.5 3.5]);
%set(gca, 'YScale', 'log', 'ylim', [1E-3 1E+3]);
title('Wait Time Transitions');
xlabel('Transitions');
ylabel('Wait Time (s)');
savefig('box_wait.fig');
print('-dbmp', 'box_wait.bmp'); 
disp('--> Boxplot created for wait time transitions');

% Create normalized current blockade (%) boxplot of all transitions.
figure(3);
boxplot(C, CG); 
set(gca, 'ylim', [0.0 0.7001]);
title('Normalized Current Blockade Level (%)');
xlabel('Transitions');
ylabel('Normalized Current Blockade Level (%)');
savefig('box_ncurr.fig');
print('-dbmp', 'box_ncurr.bmp'); 
disp('--> Boxplot created for normalized current transitions');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        BOXPLOT ANALYSIS FINALIZATION                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> Boxplots for tag transition analysis end');
fprintf('\n');
