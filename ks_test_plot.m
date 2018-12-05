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

function ks_test_plot

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
%                     CALCULATE P-VALUES USING KS TEST                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dwell time for tag capture band only (DEB/DBB).
disp('--> DWELL TIME P-VALUES USING KS TEST: TAG CAPTURE BAND ONLY (DEB/DBB)');

% Figure counter.
f = 0;
            
for i = 1:length(list)
    for j = 1:length(list)
   
        % Compute p-value of KS test.
        [h3(i,j), p3(i,j), ks2stat3(i,j)] = kstest2(S(i).DEB, S(j).DBB);
        [h, p] = kstest2(S(i).DEB, S(j).DBB);
        
        if i == j
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                                                         %
            %                           HISTOGRAM PLOTTING                            %
            %                                                                         %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                              
            disp(['--> Matrix index: ', num2str(i), '-', num2str(j)]);
        
            % Normalize histogram counts to 1 (background).
            if i == 1
                scalar = 4.8;
            elseif i == 2
                scalar = 8;
            elseif i == 3
                scalar = 2.4;
            else
                scalar = 18.75;
            end         
                        
            % Obtain number of capture events in the tag capture band only.
            [r, c] = size(S(i).DEB);

            % Display the number of maximum events captured by a single pore.
            disp(['--> Number of capture events in the tag capture band only: ', num2str(r)]);

            % Create histogram of normalization current blockade level with Gaussian fit. 
            figure(f+1);
            scalar
            XE = logspace(-3, 3, 100);
            hist([S(i).DEB; S(i).DEB; S(i).DEB; S(i).DEB], XE);
            set(gca, 'XScale', 'log');
            title('Capture events - tag capture band only');
            xlabel('Dwell time (s)');
            ylabel('Count (#)');
            xlim([1e-3 1e+3]);
            h = findobj(gca, 'Type', 'patch');
            set(h, 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'black');
            grid;
            savefig(['dwell_dist_E-', num2str(i), '.fig']);
            print('-dbmp', ['dwell_dist_E-', num2str(i), '.bmp']); 

            % Obtain number of background events in the tag capture band only.
            [r, c] = size(S(i).DBB);

            % Display the number of maximum events captured by a single pore.
            disp(['--> Number of background events in the tag capture band only: ', num2str(r)]);

            % Create histogram of normalization current blockade level with Gaussian fit. 
            figure(f+2);
            XE = logspace(-3, 3, 100);
            hist(S(i).DBB, XE);
            set(gca, 'XScale', 'log');
            title('Background events - tag capture band only');
            xlabel('Dwell time (s)');
            ylabel('Count (#)');
            xlim([1e-3 1e+3]);
            h = findobj(gca, 'Type', 'patch');
            set(h, 'FaceColor', [0.7 0.7 1], 'EdgeColor', 'black');
            grid;
            savefig(['dwell_dist_B-', num2str(i), '.fig']);
            print(['dwell_dist_B-', num2str(i), '.bmp']); 
            
            % Update figure counter.
            f = f + 2;
            
        end
                
    end
end

h3
p3
ks2stat3

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
