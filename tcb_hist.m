%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: February 16, 2016.
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: Given 4 'stats.mat' structure of static capture vs. background
% data sets, this program iterates through all of them and generates the
% dwell time histograms for:
%
% (1) DEB := dwell time of events in tag cature band (exp)
% (2) DBB := dwell time of events in tag capture band (bgr)
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function tcb_hist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         EVENT ANALYSIS STARTUP                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all')

fprintf('\n');
disp('--> Tag capture band dwell time histogram start');
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
S = struct('DEB', {}, 'DBB', {});

for i = 1:length(list)

  disp(['--> Processing structure: ', list(i).name]);  

  % Load the stats data structure into workspace.
  load(list(i).name);

  % Generate container arrays for events in tag capture band(exp).
  S(i).DEB = dwellt_tcb_exp;
  
  % Generate container arrays for events in tag capture band (bgr).
  S(i).DBB = dwellt_tcb_bgr; 

end

fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     GENERATE DWELL TIME HISTOGRAMS                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dwell time for tag capture band only (DEB/DBB).
disp('--> DWELL TIME HISTOGRAMS: TAG CAPTURE BAND ONLY (DEB/DBB)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             NUCLEOTIDE G                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> Matrix index: 1-1');

% Obtain number of capture events in the tag capture band only.
[r, c] = size(S(1).DEB);

% Display the number of maximum events captured by a single pore.
disp(['--> Number of capture events in G tag capture band only: ', num2str(r)]);

% Create histogram of normalization current blockade level with Gaussian fit. 
figure(1);
XE = logspace(-3, 3, 100);
hist([S(1).DEB; S(1).DEB; S(1).DEB; S(1).DEB; S(1).DEB], XE);
set(gca, 'XScale', 'log');
title('Capture events - G tag capture band only');
xlabel('Dwell time (s)');
ylabel('Count (#)');
xlim([1e-3 1e+3]);
h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black');
grid;
savefig('dwell_dist_E-G.fig');
print('-dbmp', 'dwell_dist_E-G.bmp'); 

% Obtain number of background events in the tag capture band only.
[r, c] = size(S(1).DBB);

% Display the number of maximum events captured by a single pore.
disp(['--> Number of background events in G tag capture band only: ', num2str(r)]);

% Create histogram of normalization current blockade level with Gaussian fit. 
figure(2);
XE = logspace(-3, 3, 100);
hist(S(1).DBB, XE);
set(gca, 'XScale', 'log');
title('Background events - G tag capture band only');
xlabel('Dwell time (s)');
ylabel('Count (#)');
xlim([1e-3 1e+3]);
h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [0.7 0.7 1], 'EdgeColor', 'black');
grid;
savefig('dwell_dist_B-G.fig');
print('dwell_dist_B-G.bmp'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             NUCLEOTIDE A                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> Matrix index: 2-2');      

% Obtain number of capture events in the tag capture band only.
[r, c] = size(S(2).DEB);

% Display the number of maximum events captured by a single pore.
disp(['--> Number of capture events in A tag capture band only: ', num2str(r)]);

% Create histogram of normalization current blockade level with Gaussian fit. 
figure(3);
XE = logspace(-3, 3, 100);
hist([S(2).DEB; S(2).DEB; S(2).DEB; S(2).DEB; S(2).DEB; S(2).DEB; S(2).DEB; S(2).DEB], XE);
set(gca, 'XScale', 'log');
title('Capture events - A tag capture band only');
xlabel('Dwell time (s)');
ylabel('Count (#)');
xlim([1e-3 1e+3]);
h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black');
grid;
savefig('dwell_dist_E-A.fig');
print('-dbmp', 'dwell_dist_E-A.bmp'); 

% Obtain number of background events in the tag capture band only.
[r, c] = size(S(2).DBB);

% Display the number of maximum events captured by a single pore.
disp(['--> Number of background events in A tag capture band only: ', num2str(r)]);

% Create histogram of normalization current blockade level with Gaussian fit. 
figure(4);
XE = logspace(-3, 3, 100);
hist(S(2).DBB, XE);
set(gca, 'XScale', 'log');
title('Background events - A tag capture band only');
xlabel('Dwell time (s)');
ylabel('Count (#)');
xlim([1e-3 1e+3]);
h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [0.7 0.7 1], 'EdgeColor', 'black');
grid;
savefig('dwell_dist_B-A.fig');
print('dwell_dist_B-A.bmp'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             NUCLEOTIDE C                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> Matrix index: 3-3');      

% Obtain number of capture events in the tag capture band only.
[r, c] = size(S(3).DEB);

% Display the number of maximum events captured by a single pore.
disp(['--> Number of capture events in C tag capture band only: ', num2str(r)]);

% Create histogram of normalization current blockade level with Gaussian fit. 
figure(5);
XE = logspace(-3, 3, 100);
hist([S(3).DEB; S(3).DEB], XE);
  
set(gca, 'XScale', 'log');
title('Capture events - C tag capture band only');
xlabel('Dwell time (s)');
ylabel('Count (#)');
xlim([1e-3 1e+3]);
h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black');
grid;
savefig('dwell_dist_E-C.fig');
print('-dbmp', 'dwell_dist_E-C.bmp'); 

% Obtain number of background events in the tag capture band only.
[r, c] = size(S(3).DBB);

% Display the number of maximum events captured by a single pore.
disp(['--> Number of background events in C tag capture band only: ', num2str(r)]);

% Create histogram of normalization current blockade level with Gaussian fit. 
figure(6);
XE = logspace(-3, 3, 100);
hist(S(3).DBB, XE);
set(gca, 'XScale', 'log');
title('Background events - C tag capture band only');
xlabel('Dwell time (s)');
ylabel('Count (#)');
xlim([1e-3 1e+3]);
h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [0.7 0.7 1], 'EdgeColor', 'black');
grid;
savefig('dwell_dist_B-C.fig');
print('dwell_dist_B-C.bmp'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             NUCLEOTIDE T                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> Matrix index: 4-4');      

% Obtain number of capture events in the tag capture band only.
[r, c] = size(S(4).DEB);

% Display the number of maximum events captured by a single pore.
disp(['--> Number of capture events in T tag capture band only: ', num2str(r)]);

% Create histogram of normalization current blockade level with Gaussian fit. 
figure(7);
XE = logspace(-3, 3, 100);
hist([S(4).DEB; S(4).DEB; S(4).DEB; S(4).DEB; S(4).DEB; S(4).DEB; ...
      S(4).DEB; S(4).DEB; S(4).DEB; S(4).DEB; S(4).DEB; S(4).DEB; ...
      S(4).DEB; S(4).DEB; S(4).DEB; S(4).DEB; S(4).DEB; S(4).DEB; ...
      S(4).DEB], XE);
  
set(gca, 'XScale', 'log');
title('Capture events - T tag capture band only');
xlabel('Dwell time (s)');
ylabel('Count (#)');
xlim([1e-3 1e+3]);
h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black');
grid;
savefig('dwell_dist_E-T.fig');
print('-dbmp', 'dwell_dist_E-T.bmp'); 

% Obtain number of background events in the tag capture band only.
[r, c] = size(S(4).DBB);

% Display the number of maximum events captured by a single pore.
disp(['--> Number of background events in T tag capture band only: ', num2str(r)]);

% Create histogram of normalization current blockade level with Gaussian fit. 
figure(8);
XE = logspace(-3, 3, 100);
hist(S(4).DBB, XE);
set(gca, 'XScale', 'log');
title('Background events - T tag capture band only');
xlabel('Dwell time (s)');
ylabel('Count (#)');
xlim([1e-3 1e+3]);
h = findobj(gca, 'Type', 'patch');
set(h, 'FaceColor', [0.7 0.7 1], 'EdgeColor', 'black');
grid;
savefig('dwell_dist_B-T.fig');
print('dwell_dist_B-T.bmp'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       TAG CAPTURE BAND FINALIZATION                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> Tag capture band dwell time histogram end');
fprintf('\n');
