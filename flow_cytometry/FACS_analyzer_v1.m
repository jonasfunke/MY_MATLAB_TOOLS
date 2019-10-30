%%
close all, clear all, clc

[filenames, pathname]=uigetfile('*.fcs','Select the fcs files','MultiSelect','on');

prefix_out = [ datestr(now, 'yyyy-mm-dd_HH-MM') '_analysis'];
path_out = [pathname prefix_out filesep];
mkdir(path_out);

%% make a cell array even if you have one file only
if ~iscell(filenames)
    filenames = {filenames};
end
   
data(1) = load_fcs_data(pathname, filenames{1}, path_out);

for i=2:length(filenames)
    data(i) = load_fcs_data(pathname, filenames{i}, path_out, data(1).roi_position);
end

%%
i_fl_ch=14;
sample_names = cell(length(filenames),1);

for j=1:length(filenames)
        sample_names{j} = [ filenames{j}(12:end-4) ' ' data(j).fcshdr.par(i_fl_ch).name];
end

%% Get limits
i = 14; % FL5-A

tmp = [];
for j=1:length(filenames)
    tmp = [tmp; data(j).fcsdat(data(j).i_gated,i_fl_ch)];
end
xlim = [min(tmp) max(tmp)];
xlim2 = [min(tmp)  median(tmp)+3*std(tmp)];
xlim3 = [min(tmp) prctile(tmp, 98)];

%% PLot FL5-A channel histogram, log scale
cur_fig = figure(1); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 20 10 ], 'PaperSize', [20 10] );

legend_tmp = {};
for j=1:length(filenames)
    tmp = data(j).fcsdat(data(j).i_gated,i);
    
   
    x=logspace(-1,8,100); % create bin edges with logarithmic scale
    histogram(tmp, x), hold on %, 'Normalization', 'pdf'
    
    legend_tmp = [legend_tmp; {[ filenames{j}(12:end-4) ' ' data(j).fcshdr.par(i_fl_ch).name]}];
    disp([filenames{j} ', ' num2str(median(tmp))])
    
end
legend(legend_tmp, 'Location', 'best')
set(gca, 'xscale','log', 'xlim', [0 xlim(2)])
xlabel(data(1).fcshdr.par(i).name), ylabel('Number of cells')

print(cur_fig, '-dpdf', [path_out filesep prefix_out '_histogram_log_FL.pdf']); %save figure

%% PLot FL5-A channel histogram

cur_fig = figure(2); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 20 10 ], 'PaperSize', [20 10] );

x_start= xlim3(1);
x_stop = xlim3(2);
h=100;
dx=20;

legend_tmp = {};
for j=1:length(filenames)
    tmp = data(j).fcsdat(data(j).i_gated,i_fl_ch);
    
   
    %x=logspace(-1,5,100); % create bin edges with logarithmic scale
    %histogram(tmp, 'Normalization', 'pdf'); hold on %, 'Normalization', 'pdf'
    [n, p, x_points] = uniform_kernel_density( tmp, h, x_start, x_stop, dx);
    plot(x_points, p), hold on
    
    legend_tmp = [legend_tmp; {[ filenames{j}(12:end-4) ' ' data(j).fcshdr.par(i_fl_ch).name]}];
    disp([filenames{j} ', ' num2str(median(tmp))])
    set(gca,  'xlim', xlim3)
    pause
end
legend(legend_tmp, 'Location', 'best')
xlabel(data(1).fcshdr.par(i_fl_ch).name), ylabel('Probability density')

print(cur_fig, '-dpdf', [path_out filesep prefix_out '_histogram_lin.pdf']); %save figure

%% plot bar graph

cur_fig = figure(3); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 10 10 ], 'PaperSize', [10 10] );

N_sample = length(filenames);
val_median = zeros(N_sample,1);
val_mean = zeros(N_sample,1);
val_err = zeros(N_sample,1);

for j=1:length(filenames)
    val_median(j) = median(data(j).fcsdat(data(j).i_gated,i_fl_ch));    
    val_mean(j) = mean(data(j).fcsdat(data(j).i_gated,i_fl_ch));    
    val_err(j) = 3*std(data(j).fcsdat(data(j).i_gated,i_fl_ch))/sqrt(length(data(j).fcsdat(data(j).i_gated,i_fl_ch)));    
    
end



bar(1:N_sample,val_median), hold on
set(gca, 'Xtick', 1:N_sample, 'XTickLabel', sample_names)
xtickangle(45)
ylabel('Median Fluorescence')
%set(gca,  'ylim', [0 1.5e4 ])%1.1*max(val_median)])
set(gca,  'ylim', [0 1.1*max(val_median)])
grid on

print(cur_fig, '-dpdf', [path_out filesep prefix_out '_median2.pdf']); %save figure



%% export csv
file_out = [path_out prefix_out '_data.txt'];
fileID = fopen(file_out,'w');

fprintf(fileID,'Name\t');
fprintf(fileID,'Median\t');
fprintf(fileID,'Mean\t');
fprintf(fileID,'\n');

for j=1:length(filenames)
    fprintf(fileID,'%s\t',sample_names{j});
    fprintf(fileID,'%f\t', median(data(j).fcsdat(data(j).i_gated,i_fl_ch)));
    fprintf(fileID,'%f\t', mean(data(j).fcsdat(data(j).i_gated,i_fl_ch)));
    fprintf(fileID,'\n');
end
fclose(fileID);

%% save data
save([path_out prefix_out '_data.mat'])

%%
t= [0.5 1 3]; %h

tmp = parula();
cc = tmp(32:64:end,:);%varycolor(3);

cur_fig = figure(4); clf

j_color = 1;
myleg = cell(0,1);
for i=[1 4 7 10]
    plot(t, val_median(i:i+2), '.-', 'Color', cc(j_color,:)), hold on
    j_color = j_color +1;
    myleg = [myleg, filenames{i}(12:end-7) ]
end

j_color = 1;
for i=[13 16 19 22]
    plot(t, val_median(i:i+2), '.--', 'Color', cc(j_color,:)), hold on
    j_color = j_color +1;
    myleg = [myleg, filenames{i}(12:end-7) ]
end
grid on
set(gca, 'YLim', [0 1.1*max(val_median)], 'XLim', [0 t(end)+1])
legend(myleg, 'location', 'best')
xlabel('Time (h)'), ylabel('Median Fluorescence')
    

set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 30 20 ], 'PaperSize', [30 20] );
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_vs_t.pdf']); %save figure


%% Plot histogram for all channels


cur_fig = figure(5); clf
N_channel = size(data(j).fcsdat,2)/2-1;

legend_tmp = {};

for j=1:length(filenames)
    
    for i=2:2:size(data(j).fcsdat,2)-2
        
        subplot(N_channel, 1, i/2)
        tmp = data(j).fcsdat(data(j).i_gated,i);
   
        x=logspace(0,7,100); % create bin edges with logarithmic scale
        histogram(tmp, x), hold on %, 'Normalization', 'pdf'
        title([ filenames{j}(12:end-4) ' ' data(j).fcshdr.par(i).name])
        %disp([filenames{j} ', ' num2str(median(tmp))])
        set(gca, 'xscale','log')
        grid on

    end
    legend_tmp = [legend_tmp; {filenames{j}(12:end-4) }];

    legend(legend_tmp, 'location', 'best')

end
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 20 10*N_channel ], 'PaperSize', [20 10*N_channel] );
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_histogram_log.pdf']); %save figure

%% 3D scatter
i=2; j=4; k=14; %FSC-A vs SSC-A vs. FL5-A

cur_fig = figure(6); clf

for l=1:3:length(filenames)
    scatter3(data(l).fcsdat(:,i), data(l).fcsdat(:,j), data(l).fcsdat(:,k), 1, '.'), hold on
    scatter3(data(l+1).fcsdat(:,i), data(l+1).fcsdat(:,j), data(l+1).fcsdat(:,k), 1, '.')
    scatter3(data(l+2).fcsdat(:,i), data(l+2).fcsdat(:,j), data(l+2).fcsdat(:,k), 1, '.')
    set(gca, 'xscale','log', 'yscale','log')
    set(gca, 'zscale','log')
    xlabel(data(l).fcshdr.par(i).name)
    ylabel(data(l).fcshdr.par(j).name)
    zlabel(data(l).fcshdr.par(k).name)
    %legend({filenames{l}(12:end-4)  filenames{l+1}(12:end-4) })
    pause
    hold off
end

%%


cur_fig = figure(6); clf

for l=1:2:length(filenames)
    scatter3(data(l).fcsdat(data(l).i_gated,i), data(l).fcsdat(data(l).i_gated,j), data(l).fcsdat(data(l).i_gated,k), 1, '.'), hold on
    scatter3(data(l+1).fcsdat(data(l+1).i_gated,i), data(l+1).fcsdat(data(l+1).i_gated,j), data(l+1).fcsdat(data(l+1).i_gated,k), 1, '.')
    set(gca, 'xscale','log', 'yscale','log')
    set(gca, 'zscale','log')
    xlabel(data(l).fcshdr.par(i).name)
    ylabel(data(l).fcshdr.par(j).name)
    zlabel(data(l).fcshdr.par(k).name)
    legend({filenames{l}(12:end-4)  filenames{l+1}(12:end-4) })
    pause
    hold off
end
