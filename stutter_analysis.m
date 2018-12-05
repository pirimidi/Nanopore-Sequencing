%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: September 9, 2015.
%
% For: Single molecule DNA sequencing via aHL nanopore array at the Church
% Lab - Genetics Department, Harvard Medical School.
%
% Purpose: Given a 'var.mat' structure of a 'vardat' run for N data sets,
% this program iterates through all of them and extracts filter level 2 
% (all clean pores) event information of current blockade level and 
% corresponding dwell time pairs, then generates a transition map of:
%
% (1) N -> N (4 distincs cases)
% (2) N -> N+1 (12 distinct cases) 
%
% calculating the frequency of each as well as the lumped (1) and (2)
% transitions.
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

function [S] = stutter_analysis(filter_e, order_e, tag_capture, dwell_cutoff)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         STUTTER ANALYSIS STARTUP                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all')

fprintf('\n');
disp('--> Stutter analysis start');
fprintf('\n');

% Set default number formatting.
format short;

% Define current working directory.
work_dir = pwd;

% Define direcory to hold figures.
if ~exist('plots', 'dir')
  mkdir('plots');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  'EXPERIMENT' STUTTER ANALYSIS SECTION                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> EXPERIMENT STUTTER STATISTICS SECTION');

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

% Define array container for wait time data.
WE = [];

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
  wait = cap.wait(:, :);  
  imed = cap.imed(:, :);

  % Filter out NaNs and diplay 2D scatter plot and histogram with a 
  % distribution fit.
  iD = ~isnan(dwell);
  
  de = dwell(iD(:));
  we = wait(iD(:));
  ce = imed(iD(:));

  % Obtain size of maximum events.
  [r, c] = size(de);

  % Display the number of maximum events captured by a single pore.
  disp(['--> Maximum number of events: ', num2str(r)]);
  fprintf('\n');

  % Collect all filtered events and update storage array.
  DE = [DE; de];
  WE = [WE; we];
  CE = [CE; ce];

end

% Obtain size of all filtered events.
[r, c] = size(DE);

% Display the number of maximum events captured by a single pore.
disp(['--> Number of all filtered events: ', num2str(r)]);
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   TAG CAPTURE STUTTER ANALYSIS SECTION                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--> TAG CAPTURE STUTTER ANALYSIS SECTION');

% Define array container for state transitions of 'experiment' data sets = 
% {G, A, C, T}.
D_GA = []; D_GC = []; D_GT = []; D_GG = [];
W_GA = []; W_GC = []; W_GT = []; W_GG = [];
C_GA = []; C_GC = []; C_GT = []; C_GG = [];

D_AA = []; D_AC = []; D_AT = []; D_AG = [];
W_AA = []; W_AC = []; W_AT = []; W_AG = [];
C_AA = []; C_AC = []; C_AT = []; C_AG = [];

D_CA = []; D_CC = []; D_CT = []; D_CG = [];
W_CA = []; W_CC = []; W_CT = []; W_CG = [];
C_CA = []; C_CC = []; C_CT = []; C_CG = [];

D_TA = []; D_TC = []; D_TT = []; D_TG = [];
W_TA = []; W_TC = []; W_TT = []; W_TG = [];
C_TA = []; C_TC = []; C_TT = []; C_TG = [];

% Define array container for all state transitions.
D_TR = []; C_TR = []; W_TR = [];

% Counter of events in tag capture box.
u = 0;

% Counter of G/A/C/T-tag capture event transitions.
ga = 1; gc = 1; gt = 1; gg = 1;
aa = 1; ac = 1; at = 1; ag = 1;
ca = 1; cc = 1; ct = 1; cg = 1;
ta = 1; tc = 1; tt = 1; tg = 1;

% Filter out events in the 'tag capture' for 'experiment' data set.
for m = 1:length(DE)

    % If event is in dwell time cutoff range (x) as well as in capture 
    % range (y), i.e., in tag capture box. 
    if (dwell_cutoff(1) < DE(m) && DE(m) < dwell_cutoff(2)) && ...
       (tag_capture(1) < CE(m) && CE(m) < tag_capture(8))  
     
        % Increment events in tag capture box.
        u = u+1;
        
        % Remember current state.
        D_TR(u) = DE(m);
        W_TR(u) = WE(m);
        C_TR(u) = CE(m);
              
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                         N->G TRANSITIONS                        %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % If current event is in G-tag capture range.
        if tag_capture(1) < C_TR(u) && C_TR(u) < tag_capture(2) 
            
            % If not the first tag capture event, i.e., will have a
            % previous event.
            if u ~= 1
                
                % If previous event is in G-tag capture range.
                if tag_capture(1) < C_TR(u-1) && C_TR(u-1) < tag_capture(2) 

                    % Remember G->G transitions.
                    D_GG(gg) = D_TR(u);
                    W_GG(gg) = W_TR(u);
                    C_GG(gg) = C_TR(u);

                    % Increment G->G transitions.
                    gg = gg+1;
                    
                end      
                
                % If previous event is in A-tag capture range.
                if tag_capture(3) < C_TR(u-1) && C_TR(u-1) < tag_capture(4) 

                    % Remember A->G transitions.
                    D_AG(ag) = D_TR(u);
                    W_AG(ag) = W_TR(u);
                    C_AG(ag) = C_TR(u);

                    % Increment A->G transitions.
                    ag = ag+1;
                    
                end    
                
                % If previous event is in C-tag capture range.
                if tag_capture(5) < C_TR(u-1) && C_TR(u-1) < tag_capture(6) 

                    % Remember C->G transitions.
                    D_CG(cg) = D_TR(u);
                    W_CG(cg) = W_TR(u);
                    C_CG(cg) = C_TR(u);

                    % Increment C->G transitions.
                    cg = cg+1;
                    
                end                 

                % If previous event is in T-tag capture range.
                if tag_capture(7) < C_TR(u-1) && C_TR(u-1) < tag_capture(8) 

                    % Remember T->G transitions.
                    D_TG(tg) = D_TR(u);
                    W_TG(tg) = W_TR(u);
                    C_TG(tg) = C_TR(u);

                    % Increment T->G transitions.
                    tg = tg+1;
                    
                end
                
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                         N->A TRANSITIONS                        %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % If current event is in A-tag capture range.
        if tag_capture(3) < CE(m) && CE(m) < tag_capture(4) 
            
            % If not the first tag capture event, i.e., will have a
            % previous event.
            if u ~= 1
                
                % If previous event is in G-tag capture range.
                if tag_capture(1) < C_TR(u-1) && C_TR(u-1) < tag_capture(2) 

                    % Remember G->A transitions.
                    D_GA(ga) = D_TR(u);
                    W_GA(ga) = W_TR(u);
                    C_GA(ga) = C_TR(u);

                    % Increment G->A transitions.
                    ga = ga+1;
                    
                end      
                
                % If previous event is in A-tag capture range.
                if tag_capture(3) < C_TR(u-1) && C_TR(u-1) < tag_capture(4) 

                    % Remember A->A transitions.
                    D_AA(aa) = D_TR(u);
                    W_AA(aa) = W_TR(u);
                    C_AA(aa) = C_TR(u);

                    % Increment A->A transitions.
                    aa = aa+1;
                    
                end    
                
                % If previous event is in C-tag capture range.
                if tag_capture(5) < C_TR(u-1) && C_TR(u-1) < tag_capture(6) 

                    % Remember C->G transitions.
                    D_CA(ca) = D_TR(u);
                    W_CA(ca) = W_TR(u);
                    C_CA(ca) = C_TR(u);

                    % Increment C->G transitions.
                    ca = ca+1;
                    
                end                 

                % If previous event is in T-tag capture range.
                if tag_capture(7) < C_TR(u-1) && C_TR(u-1) < tag_capture(8) 

                    % Remember T->A transitions.
                    D_TA(ta) = D_TR(u);
                    W_TA(ta) = W_TR(u);
                    C_TA(ta) = C_TR(u);

                    % Increment T->A transitions.
                    ta = ta+1;
                    
                end
                
            end           
        end
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                         N->C TRANSITIONS                        %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % If current event is in C-tag capture range.
        if tag_capture(5) < CE(m) && CE(m) < tag_capture(6) 
            
            % If not the first tag capture event, i.e., will have a
            % previous event.
            if u ~= 1
                
                % If previous event is in G-tag capture range.
                if tag_capture(1) < C_TR(u-1) && C_TR(u-1) < tag_capture(2) 

                    % Remember G->C transitions.
                    D_GC(gc) = D_TR(u);
                    W_GC(gc) = W_TR(u);
                    C_GC(gc) = C_TR(u);

                    % Increment G->C transitions.
                    gc = gc+1;
                    
                end      
                
                % If previous event is in A-tag capture range.
                if tag_capture(3) < C_TR(u-1) && C_TR(u-1) < tag_capture(4) 

                    % Remember A->C transitions.
                    D_AC(ac) = D_TR(u);
                    W_AC(ac) = W_TR(u);
                    C_AC(ac) = C_TR(u);

                    % Increment A->C transitions.
                    ac = ac+1;
                    
                end    
                
                % If previous event is in C-tag capture range.
                if tag_capture(5) < C_TR(u-1) && C_TR(u-1) < tag_capture(6) 

                    % Remember C->C transitions.
                    D_CC(cc) = D_TR(u);
                    W_CC(cc) = W_TR(u);
                    C_CC(cc) = C_TR(u);

                    % Increment C->C transitions.
                    cc = cc+1;
                    
                end                 

                % If previous event is in T-tag capture range.
                if tag_capture(7) < C_TR(u-1) && C_TR(u-1) < tag_capture(8) 

                    % Remember T->C transitions.
                    D_TC(tc) = D_TR(u);
                    W_TC(tc) = W_TR(u);
                    C_TC(tc) = C_TR(u);

                    % Increment T->C transitions.
                    tc = tc+1;
                    
                end
                
            end 
        end
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                 %
        %                         N->T TRANSITIONS                        %
        %                                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % If current event is in T-tag capture range.
        if tag_capture(7) < CE(m) && CE(m) < tag_capture(8) 
            
            % If not the first tag capture event, i.e., will have a
            % previous event.
            if u ~= 1
                
                % If previous event is in G-tag capture range.
                if tag_capture(1) < C_TR(u-1) && C_TR(u-1) < tag_capture(2) 

                    % Remember G->T transitions.
                    D_GT(gt) = D_TR(u);
                    W_GT(gt) = W_TR(u);
                    C_GT(gt) = C_TR(u);

                    % Increment G->T transitions.
                    gt = gt+1;
                    
                end      
                
                % If previous event is in A-tag capture range.
                if tag_capture(3) < C_TR(u-1) && C_TR(u-1) < tag_capture(4) 

                    % Remember A->T transitions.
                    D_AT(at) = D_TR(u);
                    W_AT(at) = W_TR(u);
                    C_AT(at) = C_TR(u);

                    % Increment A->T transitions.
                    at = at+1;
                    
                end    
                
                % If previous event is in C-tag capture range.
                if tag_capture(5) < C_TR(u-1) && C_TR(u-1) < tag_capture(6) 

                    % Remember C->T transitions.
                    D_CT(ct) = D_TR(u);
                    W_CT(ct) = W_TR(u);
                    C_CT(ct) = C_TR(u);

                    % Increment C->T transitions.
                    ct = ct+1;
                    
                end                 

                % If previous event is in T-tag capture range.
                if tag_capture(7) < C_TR(u-1) && C_TR(u-1) < tag_capture(8) 

                    % Remember T->T transitions.
                    D_TT(tt) = D_TR(u);
                    W_TT(tt) = W_TR(u);
                    C_TT(tt) = C_TR(u);

                    % Increment T->T transitions.
                    tt = tt+1;
                    
                end
                
            end
        end        
    end  
end

% Obtain size of N->N and N->N+1 state transitions.
[re, cw] = size(DE);

% G->N state transitions.
[rgg, cgg] = size(W_GG);
[rga, cga] = size(W_GA);
[rgc, cgc] = size(W_GC);
[rgt, cgt] = size(W_GT);

% G->N state transitions.
[raa, caa] = size(W_AA);
[rac, cac] = size(W_AC);
[rat, cat] = size(W_AT);
[rag, cag] = size(W_AG);

[rcc, ccc] = size(W_CC);
[rca, cca] = size(W_CA);
[rct, cct] = size(W_CT);
[rcg, ccg] = size(W_CG);

[rtt, ctt] = size(W_TT);
[rta, cta] = size(W_TA);
[rtc, ctc] = size(W_TC);
[rtg, ctg] = size(W_TG);

% Total number of N->N state transitions.
tnn = sum([cgg caa ccc ctt]);

% Total number of N->N+1 state transitions.
tnm = sum([cga cgc cgt cac cat cag cca cct ccg cta ctc ctg]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             STUTTER FREQUENCY                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--> Number of all filtered events: ', num2str(re)]);
disp(['--> Number of events in tag capture box: ', num2str(u)]);
disp(['--> Number of N->N state transitions: ', num2str(tnn)]);
disp(['--> Number of N->N+1 state transitions: ', num2str(tnm)]);
fprintf('\n');
disp(['--> Number of G->G state transitions: ', num2str(cgg)]);
disp(['--> Number of G->A state transitions: ', num2str(cga)]);
disp(['--> Number of G->C state transitions: ', num2str(cgc)]);
disp(['--> Number of G->T state transitions: ', num2str(cgt)]);
fprintf('\n');
disp(['--> Number of A->G state transitions: ', num2str(cag)]);
disp(['--> Number of A->A state transitions: ', num2str(caa)]);
disp(['--> Number of A->C state transitions: ', num2str(cac)]);
disp(['--> Number of A->T state transitions: ', num2str(cat)]);
fprintf('\n');
disp(['--> Number of C->G state transitions: ', num2str(ccg)]);
disp(['--> Number of C->A state transitions: ', num2str(cca)]);
disp(['--> Number of C->C state transitions: ', num2str(ccc)]);
disp(['--> Number of C->T state transitions: ', num2str(cct)]);
fprintf('\n');
disp(['--> Number of T->G state transitions: ', num2str(ctg)]);
disp(['--> Number of T->A state transitions: ', num2str(cta)]);
disp(['--> Number of T->C state transitions: ', num2str(ctc)]);
disp(['--> Number of T->T state transitions: ', num2str(ctt)]);
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            STUTTER PERCENTAGE                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--> Percent of N->N state transitions in tag capture box: ', num2str(tnn/(tnn+tnm)*100)]);
disp(['--> Percent of N->N+1 state transitions in tag capture box: ', num2str(tnm/(tnn+tnm)*100)]);
fprintf('\n');
disp(['--> Percent of G->G state transitions in tag capture box: ', num2str(cgg/(tnn+tnm)*100)]);
disp(['--> Percent of G->A state transitions in tag capture box: ', num2str(cga/(tnn+tnm)*100)]);
disp(['--> Percent of G->C state transitions in tag capture box: ', num2str(cgc/(tnn+tnm)*100)]);
disp(['--> Percent of G->T state transitions in tag capture box: ', num2str(cgt/(tnn+tnm)*100)]);
fprintf('\n');
disp(['--> Percent of A->G state transitions in tag capture box: ', num2str(cag/(tnn+tnm)*100)]);
disp(['--> Percent of A->A state transitions in tag capture box: ', num2str(caa/(tnn+tnm)*100)]);
disp(['--> Percent of A->C state transitions in tag capture box: ', num2str(cac/(tnn+tnm)*100)]);
disp(['--> Percent of A->T state transitions in tag capture box: ', num2str(cat/(tnn+tnm)*100)]);
fprintf('\n');
disp(['--> Percent of C->G state transitions in tag capture box: ', num2str(ccg/(tnn+tnm)*100)]);
disp(['--> Percent of C->A state transitions in tag capture box: ', num2str(cca/(tnn+tnm)*100)]);
disp(['--> Percent of C->C state transitions in tag capture box: ', num2str(ccc/(tnn+tnm)*100)]);
disp(['--> Percent of C->T state transitions in tag capture box: ', num2str(cct/(tnn+tnm)*100)]);
fprintf('\n');
disp(['--> Percent of T->G state transitions in tag capture box: ', num2str(ctg/(tnn+tnm)*100)]);
disp(['--> Percent of T->A state transitions in tag capture box: ', num2str(cta/(tnn+tnm)*100)]);
disp(['--> Percent of T->C state transitions in tag capture box: ', num2str(ctc/(tnn+tnm)*100)]);
disp(['--> Percent of T->T state transitions in tag capture box: ', num2str(ctt/(tnn+tnm)*100)]);
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       EVENT ANALYSIS FINALIZATION                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Navigate to working directory.
cd(work_dir);

% Generate MAT-file structure for containing all these statistics.
m = matfile('stutter_7.mat', 'Writable', true);

% All experiment transitions.
m.dwellt_all = DE;
m.dwellt_median_all = median(DE);
m.dwellt_mean_all = mean(DE);
m.dwellt_std_all = std(DE);

m.wait_all = WE;
m.wait_median_all = median(WE);
m.wait_mean_all = mean(WE);
m.wait_std_all = std(WE);

m.ncurr_all = CE;
m.ncurr_median_all = median(CE);
m.ncurr_mean_all = mean(CE);
m.ncurr_std_all = std(CE);

% N->N state transitions.
D_NN = [D_GG D_AA D_CC D_TT];
m.dwellt_NN = D_NN;
m.dwellt_median_NN = median(D_NN);
m.dwellt_mean_NN = mean(D_NN);
m.dwellt_std_NN = std(D_NN);

W_NN = [W_GG W_AA W_CC W_TT];
m.wait_NN = W_NN;
m.wait_median_NN = median(W_NN);
m.wait_mean_NN = mean(W_NN);
m.wait_std_NN = std(W_NN);

C_NN = [C_GG C_AA C_CC C_TT];
m.ncurr_NN = C_NN;
m.ncurr_median_NN = median(C_NN);
m.ncurr_mean_NN = mean(C_NN);
m.ncurr_std_NN = std(C_NN);

% N->N+1 state transitions.
D_NM = [D_GA D_GC D_GT D_AG D_AC D_AT D_CG D_CA D_CT D_TG D_TA D_TC];
m.dwellt_NM = D_NM;
m.dwellt_median_NM = median(D_NM);
m.dwellt_mean_NM = mean(D_NM);
m.dwellt_std_NM = std(D_NM);

W_NM = [W_GA W_GC W_GT W_AG W_AC W_AT W_CG W_CA W_CT W_TG W_TA W_TC];
m.wait_NM = W_NM;
m.wait_median_NM = median(W_NM);
m.wait_mean_NM = mean(W_NM);
m.wait_std_NM = std(W_NM);

C_NM = [C_GA C_GC C_GT C_AG C_AC C_AT C_CG C_CA C_CT C_TG C_TA C_TC];
m.ncurr_NM = C_NM;
m.ncurr_median_NM = median(C_NM);
m.ncurr_mean_NM = mean(C_NM);
m.ncurr_std_NM = std(C_NM);

% A->N state transitions.
m.dwellt_AA = D_AA; m.wait_AA = W_AA; m.ncurr_AA = C_AA;
m.dwellt_AC = D_AC; m.wait_AC = W_AC; m.ncurr_AC = C_AC;
m.dwellt_AT = D_AT; m.wait_AT = W_AT; m.ncurr_AT = C_AT;
m.dwellt_AG = D_AG; m.wait_AG = W_AG; m.ncurr_AG = C_AG;

% C->N state transitions.
m.dwellt_CA = D_CA; m.wait_CA = W_CA; m.ncurr_CA = C_CA;
m.dwellt_CC = D_CC; m.wait_CC = W_CC; m.ncurr_CC = C_CC;
m.dwellt_CT = D_CT; m.wait_CT = W_CT; m.ncurr_CT = C_CT;
m.dwellt_CG = D_CG; m.wait_CG = W_CG; m.ncurr_CG = C_CG;

% T->N state transitions.
m.dwellt_TA = D_TA; m.wait_TA = W_TA; m.ncurr_TA = C_TA;
m.dwellt_TC = D_TC; m.wait_TC = W_TC; m.ncurr_TC = C_TC;
m.dwellt_TT = D_TT; m.wait_TT = W_TT; m.ncurr_TT = C_TT;
m.dwellt_TG = D_TG; m.wait_TG = W_TG; m.ncurr_TG = C_TG;

% G->N state transitions.
m.dwellt_GA = D_GA; m.wait_GA = W_GA; m.ncurr_GA = C_GA;
m.dwellt_GC = D_GC; m.wait_GC = W_GC; m.ncurr_GC = C_GC;
m.dwellt_GT = D_GT; m.wait_GT = W_GT; m.ncurr_GT = C_GT;
m.dwellt_GG = D_GG; m.wait_GG = W_GG; m.ncurr_GG = C_GG;

% Return the percent occurence matrix for all transitions.
S = [cgg cga cgc cgt; cag caa cac cat; ccg cca ccc cct; ctg cta ctc ctt]/(tnn+tnm)*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%              GENERATE HISTOGRAMS OF TRANSITION FREQUENCY                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define x-axis resolution.
XE = logspace(-3, 3, 100);

%----------------------------- ALL STATES --------------------------------%

% Plot all state transitions - wait time.
figure(1);
histogram(WE, XE);
title('All state transition wait time frequency');
xlabel('Wait Time (s)');
ylabel('Frequency (cnt)');
grid;
xlim([1e-3 1e+3]);
set(gca, 'xscale', 'log');
savefig('all_wait.fig');
print('-dbmp', 'all_wait.bmp');

% Plot all state transitions - dwell time.
figure(2);
histogram(DE, XE);
title('All state transition dwell time frequency');
xlabel('Dwell Time (s)');
ylabel('Frequency (cnt)');
grid;
xlim([1e-3 1e+3]);
set(gca, 'xscale', 'log');
savefig('all_dwell.fig');
print('-dbmp', 'all_dwell.bmp');

% Plot all state transitions.
figure(3);
hist3([WE, DE], {XE XE});
title('All state transition frequency');
xlabel('Wait Time (s)');
ylabel('Dwell Time (s)');
xlim([1e-3 1e+3]);
ylim([1e-3 1e+3]);
set(gca, 'xscale', 'log', 'yscale', 'log');
set(get(gca,'child'), 'FaceColor', 'interp', 'CDataMode', 'auto');
view(2);
colormap('jet');
colorbar;
savefig('all_trans.fig');
print('-dbmp', 'all_trans.bmp');

%---------------------------- N->N STATES --------------------------------%

% Plot N->N state transitions - wait time.
figure(4);
histogram(W_NN, XE);
title('N->N state transition wait time frequency');
xlabel('Wait Time (s)');
ylabel('Frequency (cnt)');
grid;
xlim([1e-3 1e+3]);
set(gca,'xscale','log') 
savefig('NN_wait.fig');
print('-dbmp', 'NN_wait.bmp');

% Plot N->N state transitions - dwell time.
figure(5);
histogram(D_NN, XE);
title('N->N state transition dwell time frequency');
xlabel('Dwell Time (s)');
ylabel('Frequency (cnt)');
grid;
xlim([1e-3 1e+3]);
set(gca,'xscale','log') 
savefig('NN_dwell.fig');
print('-dbmp', 'NN_dwell.bmp');

% Plot all state transitions.
figure(6);
hist3([W_NN', D_NN'], {XE XE});
title('N->N state transition frequency');
xlabel('Wait Time (s)');
ylabel('Dwell Time (s)');
xlim([1e-3 1e+3]);
ylim([1e-3 1e+3]);
set(gca, 'xscale', 'log', 'yscale', 'log');
set(get(gca,'child'), 'FaceColor', 'interp', 'CDataMode', 'auto');
view(2);
colormap('jet');
colorbar;
savefig('NN_trans.fig');
print('-dbmp', 'NN_trans.bmp');

%---------------------------- N->N+1 STATES ------------------------------%

% Plot N->N+1 state transitions - wait time.
figure(7);
histogram(W_NM, XE);
title('N->N+1 state transition wait time frequency');
xlabel('Wait Time (s)');
ylabel('Frequency (cnt)');
grid;
xlim([1e-3 1e+3]);
set(gca,'xscale','log') 
savefig('NM_wait.fig');
print('-dbmp', 'NM_wait.bmp');

% Plot N->N+1 state transitions - dwell time.
figure(8);
histogram([D_GA D_GC D_GT D_AG D_AC D_AT D_CG D_CA D_CT D_TG D_TA D_TC], XE);
title('N->N+1 state transition dwell time frequency');
xlabel('Dwell Time (s)');
ylabel('Frequency (cnt)');
grid;
xlim([1e-3 1e+3]);
set(gca,'xscale','log') 
savefig('NM_dwell.fig');
print('-dbmp', 'NM_dwell.bmp');

% Plot all state transitions.
figure(9);
hist3([W_NM', ...
       [D_GA D_GC D_GT D_AG D_AC D_AT D_CG D_CA D_CT D_TG D_TA D_TC]'], {XE XE});
   
title('N->N+1 state transition frequency');
xlabel('Wait Time (s)');
ylabel('Dwell Time (s)');
xlim([1e-3 1e+3]);
ylim([1e-3 1e+3]);
set(gca, 'xscale', 'log', 'yscale', 'log');
set(get(gca,'child'), 'FaceColor', 'interp', 'CDataMode', 'auto');
view(2);
colormap('jet');
colorbar;
savefig('NM_trans.fig');
print('-dbmp', 'NM_trans.bmp');

%------------------------------ 4X4 STATES -------------------------------%
 
% Base alphabet set.
B = ['G' 'A' 'C' 'T'];

% Wait time argument set.
W = {W_GG W_GA W_GC W_GT W_AG W_AA W_AC W_AT ...
     W_CG W_CA W_CC W_CT W_TG W_TA W_TC W_TT};

% Dwell time argument set.
D = {D_GG D_GA D_GC D_GT D_AG D_AA D_AC D_AT ...
     D_CG D_CA D_CC D_CT D_TG D_TA D_TC D_TT};
 
% State iteration counter.
k = 1;

% Figure iteration counter.
f = 1;

for p = 1:length(B)
    for q = 1:length(B)
               
        % Check if the event counter is non-empty.
        if isempty(W{k}) == 0
        
            % Plot all N->M state transitions for wait time, where N = M = {G, A, C, T}.
            figure(10+3*(f-1));
            histogram(W{k}, XE);
            title([num2str(B(p)) '->' num2str(B(q)) ' state transition wait time frequency']);
            xlabel('Wait Time (s)');
            ylabel('Frequency (cnt)');
            grid;
            xlim([1e-3 1e+3]);
            set(gca,'xscale','log') 
            savefig([num2str(B(p)) num2str(B(q)) '_wait.fig']);
            print('-dbmp', [num2str(B(p)) num2str(B(q)) '_wait.bmp']);     

            % Plot all N->M state transitions for dwell time, where N = M = {G, A, C, T}.
            figure(11+3*(f-1));
            histogram(D{k}, XE);
            title([num2str(B(p)) '->' num2str(B(q)) ' state transition dwell time frequency']);
            xlabel('Dwell Time (s)');
            ylabel('Frequency (cnt)');
            grid;
            xlim([1e-3 1e+3]);
            set(gca,'xscale','log') 
            savefig([num2str(B(p)) num2str(B(q)) '_dwell.fig']);
            print('-dbmp', [num2str(B(p)) num2str(B(q)) '_dwell.bmp']);

            % Plot all state transitions.
            figure(12+3*(f-1));
            hist3([W{k}', D{k}'], {XE XE});
            title([num2str(B(p)) '->' num2str(B(q)) ' state transition frequency']);
            xlabel('Wait Time (s)');
            ylabel('Dwell Time (s)');
            xlim([1e-3 1e+3]);
            ylim([1e-3 1e+3]);
            set(gca, 'xscale', 'log', 'yscale', 'log');
            set(get(gca,'child'), 'FaceColor', 'interp', 'CDataMode', 'auto');
            view(2);
            colormap('jet');
            colorbar;
            savefig([num2str(B(p)) num2str(B(q)) '_trans.fig']);
            print('-dbmp', [num2str(B(p)) num2str(B(q)) '_trans.bmp']);
            
            % Increment figure interation counter.
            f = f + 1;
        
        end 
        
        % Increment state interation counter.
        k = k + 1;
        
    end
end

% Move all figures to 'plots' directory.
movefile('*.fig', 'plots');
movefile('*.bmp', 'plots');

% Navigate to working directory.
cd(work_dir);

disp('--> Stutter analysis end');
fprintf('\n');
