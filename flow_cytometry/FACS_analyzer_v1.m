%% load fcs data
close all, clear all, clc

[filenames, pathname]=uigetfile('*.fcs','Select the fcs files','MultiSelect','on');

prefix_out = [ datestr(now, 'yyyy-mm-dd_HH-MM') '_analysis'];
path_out = [pathname prefix_out filesep];
mkdir(path_out);

%% make a cell array even if you have one file only
if ~iscell(filenames)
    filenames = {filenames};
end
   
data(1) = load_fcs_data(pathname, filenames{1}, path_out, 0.05);

for i=2:length(filenames)
    data(i) = load_fcs_data(pathname, filenames{i}, path_out, 0.05, data(1).roi_position);
    %data(i) = load_fcs_data(pathname, filenames{i}, path_out, 0.1);
end

%% create sample names
i_fl_ch=14; % FL5-A
i_fsc_ch=2; 
i_ssc_ch=4; %FSC-A vs SSC-A


sample_names = cell(length(filenames),1);

for j=1:length(filenames)
        sample_names{j} = [ filenames{j}(12:end-4) ' ' data(j).fcshdr.par(i_fl_ch).name];
end

%% Get limits for plots
tmp = [];
for j=1:length(filenames)
    tmp = [tmp; data(j).fcsdat(data(j).i_gated,i_fl_ch)];
end
xlim = [min(tmp) max(tmp)];
xlim2 = [min(tmp)  median(tmp)+3*std(tmp)];
xlim3 = [min(tmp) prctile(tmp, 98)];

%% make color scheme
tmp = parula;
cc = tmp(1:floor(size(tmp,1)/(length(filenames))):end,:);

%% PLot FL5-A channel histogram, log scale
cur_fig = figure(1); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 20 10 ], 'PaperSize', [20 10] );

legend_tmp = {};
for j=1:length(filenames)
    tmp = data(j).fcsdat(data(j).i_gated,i_fl_ch);
    
    x=logspace(-1,8,200); % create bin edges with logarithmic scale
    n = histogram(tmp, x, 'DisplayStyle','stairs', 'Normalization', 'probability', 'EdgeColor', cc(j,:)); hold on %, 'Normalization', 'pdf'
    %plot(x(1:end-1), n.Values)
    legend_tmp = [legend_tmp; {[ filenames{j}(12:end-4) ' ' data(j).fcshdr.par(i_fl_ch).name]}];
    %disp([filenames{j} ', ' num2str(median(tmp))])
    
end
legend(legend_tmp, 'Location', 'best')
set(gca, 'xscale','log', 'xlim', [0 xlim(2)])
xlabel(data(1).fcshdr.par(i_fl_ch).name), ylabel('Fraction of cells')

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
    plot(x_points, p, 'Color', cc(j,:)), hold on
    
    legend_tmp = [legend_tmp; {[ filenames{j}(12:end-4) ' ' data(j).fcshdr.par(i_fl_ch).name]}];
    %disp([filenames{j} ', ' num2str(median(tmp))])
    set(gca,  'xlim', xlim3)
 end
legend(legend_tmp, 'Location', 'best')
xlabel(data(1).fcshdr.par(i_fl_ch).name), ylabel('Probability density')

print(cur_fig, '-dpdf', [path_out filesep prefix_out '_histogram_lin.pdf']); %save figure

%% plot bar graph
cur_fig = figure(3); clf
N_sample = length(filenames);
val_median = zeros(N_sample,1);
val_median_nongated = zeros(N_sample,1);
val_mean = zeros(N_sample,1);
val_err = zeros(N_sample,1);
val_ssc = zeros(N_sample,1);
val_fsc = zeros(N_sample,1);

for j=1:length(filenames)
    val_median(j) = median(data(j).fcsdat(data(j).i_gated,i_fl_ch));    
    val_median_nongated(j) = median(data(j).fcsdat(:,i_fl_ch));    
    val_mean(j) = mean(data(j).fcsdat(data(j).i_gated,i_fl_ch));    
    val_err(j) = 3*std(data(j).fcsdat(data(j).i_gated,i_fl_ch))/sqrt(length(data(j).fcsdat(data(j).i_gated,i_fl_ch)));    

    ssc_median(j) = median(data(j).fcsdat(data(j).i_gated,i_ssc_ch));   
    fsc_median(j) = median(data(j).fcsdat(data(j).i_gated,i_fsc_ch));   
end

subplot(3, 1, 1)
plot(1:N_sample,ssc_median, '.'), hold on
plot(1:N_sample,fsc_median, '.')
legend({'Median SSC', 'Median FSC'}, 'location', 'best')
grid on
ylabel('Median SC')
set(gca, 'Xtick', 1:N_sample, 'XTickLabel', [], 'Xlim', [0 N_sample+1])

subplot(3, 1, 2:3)
bar(1:N_sample,val_median), hold on
plot(1:N_sample,val_median_nongated, '.')
set(gca, 'Xtick', 1:N_sample, 'XTickLabel', sample_names, 'Xlim', [0 N_sample+1])
xtickangle(45)
ylabel('Median Fluorescence')
%set(gca,  'ylim', [0 1.5e4 ])%1.1*max(val_median)])
set(gca,  'ylim', [0 1.1*max(val_median)])
legend({'gated', 'non-gated'}, 'Location', 'best')
grid on

set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 max(0.8*N_sample,12) max(8*N_sample*0.05, 10) ], ...
    'PaperSize', [max(0.8*N_sample,12) max(8*N_sample*0.05,10)] );

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
disp('txt file written.')

%% save data
close all
save([path_out prefix_out '_data.mat'])
disp('Data saved.')

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
        title([ filenames{j}(12:end-4) ' ' data(j).fcshdr.par(i).name ' (' num2str(i) ')'])
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


%%
cur_fig = figure(6); clf
scatter_lim = [min(data(1).fcsdat(:,2)) max(data(j).fcsdat(:,2)) ...
    min(data(1).fcsdat(:,4)) max(data(j).fcsdat(:,4))];
NN_lim = [min(data(1).NN) max(data(1).NN)];
tmp1 = [prctile(data(j).fcsdat(:,i_fl_ch),10) prctile(data(j).fcsdat(:,i_fl_ch),95)];
% 4 = SSC-A, 2 = FSC-A
for j=2:length(filenames)
    if scatter_lim(1) > min(data(j).fcsdat(:,2))
        scatter_lim(1) = min(data(j).fcsdat(:,2));
    end
    if scatter_lim(2) < max(data(j).fcsdat(:,2))
        scatter_lim(2) = max(data(j).fcsdat(:,2));
    end
    if scatter_lim(3) > min(data(j).fcsdat(:,4))
        scatter_lim(3) = min(data(j).fcsdat(:,4));
    end
    if scatter_lim(4) < max(data(j).fcsdat(:,4))
        scatter_lim(4) = max(data(j).fcsdat(:,4));
    end
    if NN_lim(2) < max(data(j).NN)
        NN_lim(2) = max(data(j).NN);
    end
    if NN_lim(1) > min(data(j).NN)
        NN_lim(1) = min(data(j).NN);
    end
    tmp1(j,:) = [prctile(data(j).fcsdat(:,i_fl_ch),10) prctile(data(j).fcsdat(:,i_fl_ch),95)];
    
    
    
    
end
Fl_lim = [min(tmp1(:,1)) min(tmp1(:,2))];
scatter_lim(scatter_lim(:)<0)=100;


scatter_path = [path_out 'scatter_plots' filesep];
mkdir(scatter_path)
% set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
%     'PaperPosition', [0 0 15 10 ], 'PaperSize', [15 10] );
% for j=1:length(filenames)
%     scatter(data(j).fcsdat(:,2),data(j).fcsdat(:,4), 10, (data(j).fcsdat(:,i_fl_ch)), '.')
%     colorbar
%     xlabel(data(j).fcshdr.par(2).name), ylabel(data(j).fcshdr.par(4).name)
%     grid on
%     caxis([prctile(data(j).fcsdat(:,i_fl_ch),10) prctile(data(j).fcsdat(:,i_fl_ch),95)])
%     h = colorbar;
%     ylabel(h, data(j).fcshdr.par(i_fl_ch).name)
%     set(gca,'xscale','log','yscale','log', 'XLim', scatter_lim(1:2) , 'YLim', scatter_lim(3:4) )
%     title(sample_names{j})
%     print(cur_fig, '-dpdf', [scatter_path filenames{j}(1:end-4) '_SSC-FSC-FL5.pdf']); %save figure
% end
% disp('done')

%%
N_column = 5; % spalten
N_row = ceil(length(filenames)/N_column); % zeilen
cur_fig = figure(8); clf
for j=1:length(filenames)
    subplot(N_row, N_column, j)
    scatter(data(j).fcsdat(:,2),data(j).fcsdat(:,4), 5, data(j).fcsdat(:,i_fl_ch), '.')
    h = colorbar;
    ylabel(h, data(j).fcshdr.par(i_fl_ch).name)
    xlabel(data(j).fcshdr.par(2).name), ylabel(data(j).fcshdr.par(4).name)
    grid on
    caxis(Fl_lim)
    set(gca,'xscale','log','yscale','log', 'XLim', scatter_lim(1:2) , 'YLim', scatter_lim(3:4) )
    title(filenames{j}(12:end-4))
end
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*14 N_row*12 ], 'PaperSize', [N_column*14 N_row*12 ] );
print(cur_fig, '-dpdf', [scatter_path 'Overview_SSC-FSC-' data(1).fcshdr.par(i_fl_ch).name  '.pdf']); %save figure


%% scatter overview

cur_fig = figure(7); clf
for j=1:length(filenames)
    subplot(N_row, N_column, j)
    scatter(data(j).fcsdat(:,2),data(j).fcsdat(:,4), 5, (data(j).NN(:)), '.')
    h = colorbar;
    ylabel(h, 'Nearest Neighbors')
    xlabel(data(j).fcshdr.par(2).name), ylabel(data(j).fcshdr.par(4).name)
    grid on
    caxis(NN_lim)
    set(gca,'xscale','log','yscale','log', 'XLim', scatter_lim(1:2) , 'YLim', scatter_lim(3:4) )
    title(filenames{j}(12:end-4))
end
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*13 N_row*12 ], 'PaperSize', [N_column*13 N_row*12 ] );
print(cur_fig, '-dpdf', [scatter_path 'Overview_SSC-FSC-NN.pdf']); %save figure
