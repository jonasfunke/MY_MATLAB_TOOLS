close all, clear all, clc

% set fluorescence channels
i_fsc_ch=2; 
i_ssc_ch=3; %FSC-A , for Cytoflex 4
ch_R1 = 4 %11; % FL5-A
%ch_R2 = 12; % FL5-A
%ch_R3 = 13; % FL5-A
i_fl_ch = ch_R1;

% load fcs data
[filenames, pathname]=uigetfile('*.fcs','Select the fcs files','MultiSelect','on');



%% create output dir
prefix_out = [ datestr(now, 'yyyy-mm-dd_HH-MM') '_analysis'];
tmp = inputdlg({'Name of analysis (prefix):'}, 'Name of analysis (prefix):' , 1, {prefix_out} );
prefix_out = tmp{1};
path_out = [pathname prefix_out filesep];
mkdir(path_out);

%% make a cell array even if you have one file only
if ~iscell(filenames)
    filenames = {filenames};
end
   
data(1) = load_fcs_data_attune(pathname, filenames{1}, path_out, 0.1);

for i=2:length(filenames)
    data(i) = load_fcs_data_attune(pathname, filenames{i}, path_out, 0.1, data(1).roi_position);
    %data(i) = load_fcs_data(pathname, filenames{i}, path_out, 0.05);
end

%% create sample names
sample_names = cell(length(filenames),1);

for j=1:length(filenames)
        sample_names{j} = filenames{j}(12:end-4);
end


%% calculate plotting limits

NN_lim = [min(data(1).NN) max(data(1).NN)];

for j=2:length(filenames)
    if NN_lim(2) < max(data(j).NN)
        NN_lim(2) = max(data(j).NN);
    end
    if NN_lim(1) > min(data(j).NN)
        NN_lim(1) = min(data(j).NN);
    end
end
scatter_lim = [1e4 2^20 1e4 2^20];

%% make fcs-ssc scatter plots
N_column = ceil(16*sqrt(length(filenames)/16/9)); % 8; % spalten
N_row = ceil(length(filenames)/N_column); % zeilen


cur_fig = figure(1); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*13 N_row*12 ], 'PaperSize', [N_column*13 N_row*12 ] );
for j=1:length(filenames)
    subplot(N_row, N_column, j)
    scatter(data(j).fcsdat(:,i_fsc_ch),data(j).fcsdat(:,i_ssc_ch), 5, (data(j).NN(:)), '.'), hold on
    cur_pgon = polyshape(data(j).roi_position);
    plot(cur_pgon, 'FaceColor', 'none', 'EdgeColor', 'r')
    legend([num2str(round(100*sum(data(j).i_gated)/length(data(j).i_gated))) '% gated (' num2str(sum(data(j).i_gated)) ' of ' num2str(length(data(j).i_gated)) ')'],'Location', 'SouthEast')
    xlabel(data(j).fcshdr.par(i_fsc_ch).name), ylabel(data(j).fcshdr.par(i_ssc_ch).name)
    
    grid on
    caxis(NN_lim)
    set(gca,'xscale','log','yscale','log', 'XLim', scatter_lim(1:2) , 'YLim', scatter_lim(3:4) )
    title(sample_names{j})
    
    if j == length(filenames)
        h = colorbar;
        ylabel(h, 'Nearest Neighbors')
    end
end

pause(1)
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_overview_SSC-FSC-NN.pdf']); %save figure

%% calculate median values
N_channel = size(data(j).fcsdat,2);
N_sample = length(filenames);

val_median = zeros(N_sample,N_channel);
val_median_nongated = zeros(N_sample,N_channel);
val_mean = zeros(N_sample,N_channel);
val_err = zeros(N_sample,N_channel);
N_counts = zeros(N_sample,2);

for j=1:length(filenames)
    for k=1:size(data(j).fcsdat,2)
        val_median(j,k) = median(data(j).fcsdat(data(j).i_gated,k));    
        val_median_nongated(j,k) = median(data(j).fcsdat(:,k));    
        val_mean(j,k) = mean(data(j).fcsdat(data(j).i_gated,k));    
        val_err(j,k) = 3*std(data(j).fcsdat(data(j).i_gated,k))/sqrt(length(data(j).fcsdat(data(j).i_gated,k)));    
    end
    %ssc_median(j) = median(data(j).fcsdat(data(j).i_gated,i_ssc_ch));   
    %fsc_median(j) = median(data(j).fcsdat(data(j).i_gated,i_fsc_ch));   
    
    N_counts(j,1) = length(data(j).fcsdat(data(j).i_gated,i_fl_ch));
    N_counts(j,2) = length(data(j).fcsdat(:,i_fl_ch));
end

%%  and make histogram for each channel



cur_fig = figure(2); clf
subplot(2, 1, 1)
plot(1:N_sample,N_counts(:,1), '.'), hold on
plot(1:N_sample,N_counts(:,2), '.')
legend({'Events (gated)', 'Events (total)'}, 'location', 'best')
grid on
ylabel('Events')
set(gca, 'Xtick', 1:N_sample, 'XTickLabel', [], 'Xlim', [0 N_sample+1])

subplot(2, 1, 2)
plot(1:N_sample,val_median(:,i_ssc_ch), '.'), hold on
plot(1:N_sample,val_median(:,i_fsc_ch), '.'), hold on
legend({'Median SSC', 'Median FSC'}, 'location', 'best')
grid on
ylabel('Median SC')
set(gca, 'Xtick', 1:N_sample, 'XTickLabel', sample_names, 'Xlim', [0 N_sample+1])
xtickangle(30)

set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 max(0.8*N_sample,12) max(8*N_sample*0.05, 20) ], ...
    'PaperSize', [max(0.8*N_sample,12) max(8*N_sample*0.05,20)] );

print(cur_fig, '-dpdf', [path_out filesep prefix_out '_counts_sc.pdf']); %save figure


for i=3:size(data(1).fcsdat,2)
    cur_fig = figure(3); clf

    bar(1:N_sample,val_median(:,i)), hold on
    plot(1:N_sample,val_median_nongated(:,i), '.')
    set(gca, 'Xtick', 1:N_sample, 'XTickLabel', sample_names, 'Xlim', [0 N_sample+1])
    xtickangle(30)
    ylabel(data(1).fcshdr.par(i).name)
    set(gca,  'ylim', [0 1.1*max(val_median(:,i))])
    legend({'gated', 'non-gated'}, 'Location', 'best')
    grid on

    set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
        'PaperPosition', [0 0 max(1*N_sample,12) 10 ], ...
        'PaperSize', [max(1*N_sample,12) 10] );

    print(cur_fig, '-dpdf', [path_out filesep prefix_out '_median_ch-' data(1).fcshdr.par(i).name '.pdf']); %save figure
end


%% save data
close all
save([path_out prefix_out '_data.mat'])
disp('Data saved.')

%%
k = 1
l = i_fsc_ch

cur_fig = figure(4); clf
for j=1:length(filenames)
    subplot(N_row, N_column, j)
    
    scatter(data(j).fcsdat(:,k),data(j).fcsdat(:,l), 5, '.'), hold on
    
    xlabel(data(j).fcshdr.par(k).name), ylabel(data(j).fcshdr.par(l).name)
    
    
    
    grid on
    %caxis(NN_lim)
    %set(gca,'xscale','log','yscale','log', 'XLim', scatter_lim(1:2) , 'YLim', scatter_lim(3:4) )
    set(gca,'yscale','log',  'YLim', scatter_lim(1:2) )
    title(sample_names{j})

end

%%
k = 3
l = 6

cur_fig = figure(4); clf
for j=1:length(filenames)
    subplot(N_row, N_column, j)
    
    scatter(data(j).fcsdat(:,k),data(j).fcsdat(:,l), 5, '.'), hold on
    scatter(data(j).fcsdat(data(j).i_gated,k),data(j).fcsdat(data(j).i_gated,l), 5, '.'), hold on
    
    xlabel(data(j).fcshdr.par(k).name), ylabel(data(j).fcshdr.par(l).name)
    
    
    
    grid on
    %caxis(NN_lim)
    %set(gca,'xscale','log','yscale','log', 'XLim', scatter_lim(1:2) , 'YLim', scatter_lim(3:4) )
    %set(gca,'xscale','log','yscale','log', 'XLim', scatter_lim(1:2) )
    set(gca,'XLim', scatter_lim(1:2) , 'YLim', scatter_lim(3:4) )
    title(sample_names{j})

end



%%
k = 2
l = 3
m = 6

j=1
cur_fig = figure(5); clf

scatter3(data(j).fcsdat(:,k),data(j).fcsdat(:,l), data(j).fcsdat(:,m), 4, '.'), hold on
scatter3(data(j).fcsdat(data(j).i_gated,k),data(j).fcsdat(data(j).i_gated,l), data(j).fcsdat(data(j).i_gated,m), 5, '.'), hold on
%set(gca,'xscale','log','yscale','log' )

xlabel(data(j).fcshdr.par(k).name)
ylabel(data(j).fcshdr.par(l).name)
zlabel(data(j).fcshdr.par(m).name)

%%
disp('Done')