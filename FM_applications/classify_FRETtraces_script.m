%%  startup
close all, clear all, clc

% load traces
tmp1 = load('/Volumes/JonasStorage/FM_data/session2/2015-05-22_18-59_analysis/data.mat');
tmp2 = load('/Volumes/JonasStorage/FM_data/session3_200ms/2015-05-22_20-20_analysis/data.mat');
tmp3 = load('/Volumes/JonasStorage/FM_data/session4_200ms/2015-05-23_11-05_analysis/data.mat');

path_out = '/Volumes/JonasStorage/FM_data/2015-05-23_inner_analysis';
%%
traces = [tmp1.traces; tmp2.traces; tmp3.traces];

%% loop through traces and compute FRET / print traces
n_traces = size(traces,1);
index = zeros(n_traces,1);
%%
cf = figure();
i = 1062;
while i <= n_traces
    
    Idd = traces{i,1}(:,4);
    Ida = traces{i,2}(:,4);
    
    x_dd = traces{i,1}(1,2)+1;
    y_dd = traces{i,1}(1,3)+1;
    x_da = traces{i,2}(1,2)+1;
    y_da = traces{i,2}(1,3)+1;
    x_aa = traces{i,3}(1,2)+1;
    y_aa = traces{i,3}(1,3)+1;
    

    subplot(3,4, 1:4)
    plot(traces{i,1}(:,1), traces{i,1}(:,4), 'g-', traces{i,2}(:,1), traces{i,2}(:,4), 'b-',  'Linewidth', 1)
     ylabel('Integr. Intensity')
    title([num2str(i) ' of ' num2str(size(traces,1))])
    
    subplot(3,4, 5:8)
    plot(traces{i,3}(:,1), traces{i,3}(:,4), 'r',  'Linewidth', 1)
    ylabel('Integr. Intensity')
        
    subplot(3,4, 9:12)
    plot(traces{i,1}(:,1), Ida ./ (Idd+Ida), 'k', 'Linewidth', 1)

    xlabel('Framenumber'), ylabel('E_app')
    set(gca, 'YLim', [min(Ida ./ (Idd+Ida)) max(Ida ./ (Idd+Ida))])
    set(gca, 'YLim', [0.3 0.6])

    print_str = questdlg('Discard or save', 'Discard or save?', 'Discard','Save and print','stop', 'Discard');
    
    if strcmp(print_str,'Save and print')
       print(cf, '-dtiff', '-r300', [path_out filesep 'trace_' sprintf('%.03i',i) '.tif'])
       index(i) = 1;
    else
        index(i) = 0;
    end
    if strcmp(print_str,'stop')
        i = n_traces+1;
    else
            i = i+1;
    end
end
index
