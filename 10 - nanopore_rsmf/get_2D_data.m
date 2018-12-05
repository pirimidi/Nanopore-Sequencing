function get_2D_data(filename)
% To obtain XY data arrays displayed in a normalized single nanopore trace. 
% Need to remove histogram panel from the composite plot manually first to
% make it work.
%
% Example usage: [x, y] = get_2D_data('example2D.fig')

    if nargin==1
        fig1 = load(filename, '-mat');

        % Obtain XY data arrays from 2D plot.
        D = fig1.hgS_070000.children.children;
        X = D(1).properties.XData;
        Y = D(1).properties.YData;

    else
        disp('--> Usage: [x, y] = get_2D_data(example2D.fig)');
    end
    
    % Plot XY data for check.
    hold on;
    plot(X, Y, 'red');

    % Calculate total number of time steps, T.
    T = length(X)
    
    % Calculate the average normalized current [0,1] for the single trace.
    I_ave = mean(Y)
    
    % Display current value for the first and last time step.
    I_first = Y(1)
    I_last = Y(T)
    
    % Current difference summation holder.
    SUM =0;
    
    % Iterate through all enrties to calculate RMSF value.
    for i = 1:length(X)
        
        I_diff = (Y(i) - I_ave)^2;
        SUM = SUM + I_diff;
        
    end
    
    SUM
    RMSF = sqrt(1/T*SUM)        
            
end
