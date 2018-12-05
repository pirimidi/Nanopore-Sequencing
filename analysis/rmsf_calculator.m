%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: January 18, 2016.
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: Given a set of single pore current traces (.fig), this program
% calculates the root-mean-square fluctuation for each trace, then
% calculates the average value for the entire set (notind STD values). 
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function rmsf_calculator(I_open)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         STUTTER ANALYSIS STARTUP                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all')

fprintf('\n');
disp('--> RMSF calculator start');
fprintf('\n');

% Set default number formatting.
format short;

% Define current working directory.
work_dir = pwd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   INDIVIDUAL RMSF CALCULATOR SECTION                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> INDIVIDUAL RMSF CALCULATOR SECTION');

% Navigate to 'pore_traces' data directory.
if ~exist('pore_traces', 'dir')
  mkdir('pore_traces');
end

cd 'pore_traces';

% Read in all 'pore trace' text file names one-by-one.
list = dir('raw_current_*');

% Define array container for average open channel current values for each 
% single pore.
I_ave_values = []; 

% Define array container for individual RMSF value.
RMSF_values = []; 

for i = 1:length(list)

    % Load in current figure data.
    fig = load(list(i).name, '-mat');
    
    disp(['--> Processing file: ', list(i).name]); 

    % Obtain XY data arrays from current current trace.
    D = fig.hgS_070000.children.children;
    X = D(1).properties.XData;
    Y = D(1).properties.YData;

    % Plot XY data for check.
    hold on;
    plot(X, Y, 'red');

    % Calculate total number of time steps, T.
    T = length(X)
    
    % Calculate the average normalized current [0,1] for the single trace.
    I_ave = I_open * mean(Y)
    
    % Current difference summation holder.
    SUM =0;
    
    % Iterate through all enrties to calculate RMSF value.
    for j = 1:length(X)
        
        I_diff = ((I_open * Y(j)) - I_ave)^2;
        SUM = SUM + I_diff;
        
    end
    
    SUM
    RMSF = sqrt(1/T*SUM)    
    
    I_ave_values(i) = I_ave;
    RMSF_values(i) = RMSF;
            
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          RMSF ANALYSIS SUMMARY                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Navigate to working directory.
cd(work_dir);

% Display the average open channel current for each analyzed single pore trace.
I_ave_values

% Display the average open channel current.
disp(['--> Average open channel current (pA): ', num2str(mean(I_ave_values))]);

% Display the standard deviation of open channel current.
disp(['--> STD of open channel current (pA): ', num2str(std(I_ave_values))]);

% Display the individual RMSF values for each analyzed single pore trace.
fprintf('\n');
RMSF_values

% Display the average RMSF value.
disp(['--> Average RMSF value (pA): ', num2str(mean(RMSF_values))]);

% Display the standard deviation of RMSF value.
disp(['--> STD of RMSF values (pA): ', num2str(std(RMSF_values))]);

% Close all opened figures.
close all;

fprintf('\n');
disp('--> RMSF calculator end');
fprintf('\n');
