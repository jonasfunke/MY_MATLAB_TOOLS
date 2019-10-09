%%
close all, clear all, clc

[filenames, pathname]=uigetfile('*.fcs','Select the fcs files','MultiSelect','on');

prefix_out = [ datestr(now, 'yyyy-mm-dd_HH-MM') '_analysis'];
path_out = [pathname prefix_out filesep];
mkdir(path_out);


%%

data(1) = load_fcs_data(pathname, filenames{1}, path_out);

for i=2:length(filenames)
    data(i) = load_fcs_data(pathname, filenames{i}, path_out, data(1).roi_position);

end

%%
i=14;
sample_names = cell(length(filenames),1);

for j=1:length(filenames)
        sample_names{j} = [ filenames{j}(12:end-4) ' ' data(j).fcshdr.par(i).name];
end

%% Get limits
i = 14; % FL5-A

tmp = [];
for j=1:length(filenames)
    tmp = [tmp; data(j).fcsdat(data(j).i_gated,i)];
end
xlim = [min(tmp) max(tmp)];
xlim2 = [min(tmp)  median(tmp)+3*std(tmp)];
xlim3 = [min(tmp) prctile(tmp, 90)];

%% PLot FL5-A channel histogram, log scale

cur_fig = figure(1); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 20 10 ], 'PaperSize', [20 10] );

legend_tmp = {};
for j=1:length(filenames)
    tmp = data(j).fcsdat(data(j).i_gated,i);
    
   
    x=logspace(-1,5,100); % create bin edges with logarithmic scale
    histogram(tmp, x), hold on %, 'Normalization', 'pdf'
    
    legend_tmp = [legend_tmp; {[ filenames{j}(12:end-4) ' ' data(j).fcshdr.par(i).name]}];
    disp([filenames{j} ', ' num2str(median(tmp))])
end
legend(legend_tmp, 'Location', 'best')
set(gca, 'xscale','log', 'xlim', [0 xlim(2)])
xlabel(data(1).fcshdr.par(i).name), ylabel('Counts')

print(cur_fig, '-dpdf', [path_out filesep prefix_out '_histogram_log.pdf']); %save figure

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
    tmp = data(j).fcsdat(data(j).i_gated,i);
    
   
    %x=logspace(-1,5,100); % create bin edges with logarithmic scale
    %histogram(tmp, 'Normalization', 'pdf'); hold on %, 'Normalization', 'pdf'
    [n, p, x_points] = uniform_kernel_density( tmp, h, x_start, x_stop, dx);
    plot(x_points, p), hold on
    
    legend_tmp = [legend_tmp; {[ filenames{j}(12:end-4) ' ' data(j).fcshdr.par(i).name]}];
    disp([filenames{j} ', ' num2str(median(tmp))])
end
legend(legend_tmp, 'Location', 'best')
set(gca,  'xlim', xlim3)
xlabel(data(1).fcshdr.par(i).name), ylabel('Probability density')

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
    val_median(j) = median(data(j).fcsdat(data(j).i_gated,i));    
    val_mean(j) = mean(data(j).fcsdat(data(j).i_gated,i));    
    val_err(j) = 3*std(data(j).fcsdat(data(j).i_gated,i))/sqrt(length(data(j).fcsdat(data(j).i_gated,i)));    
    
end



bar(1:N_sample,val_median), hold on
set(gca, 'Xtick', 1:N_sample, 'XTickLabel', sample_names)
xtickangle(45)
ylabel('Median Fluorescence')
set(gca,  'ylim', [0 1.1*max(val_median)])
grid on

print(cur_fig, '-dpdf', [path_out filesep prefix_out '_median.pdf']); %save figure



%% export csv
file_out = [path_out prefix_out '_data.txt'];
fileID = fopen(file_out,'w');
i=14;

fprintf(fileID,'Name\t');
fprintf(fileID,'Median\t');
fprintf(fileID,'Mean\t');
fprintf(fileID,'\n');

for j=1:length(filenames)
    fprintf(fileID,'%s\t',sample_names{j});
    fprintf(fileID,'%f\t', median(data(j).fcsdat(data(j).i_gated,i)));
    fprintf(fileID,'%f\t', mean(data(j).fcsdat(data(j).i_gated,i)));
    fprintf(fileID,'\n');
end
fclose(fileID);


%%




