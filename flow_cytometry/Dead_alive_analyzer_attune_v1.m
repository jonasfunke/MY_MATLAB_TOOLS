close all, clear all, clc

% set fluorescence channels
i_fsc_ch=2; 
i_ssc_ch=3; %FSC-A , for Cytoflex 4


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

for i=1:length(filenames)
    data(i) = load_fcs_data_attune_v2(pathname, filenames{i}, 0.1);
end

%% create sample names
i=1;
while i<length(filenames{1})
    pattern = filenames{1}(1:end-i)
    if all(startsWith(filenames, pattern))
        i_found=i-1;
        i=length(filenames{1}); % stop
    end
    i = i+1;
end

sample_names = cell(length(filenames),1);

for j=1:length(filenames)
        sample_names{j} = filenames{j}(end-i_found:end-4);
        disp(sample_names{j})
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

%%

xy_combined = zeros(0, 2);
for j=1:length(filenames)
    xy = [data(j).fcsdat(:,i_fsc_ch),data(j).fcsdat(:,i_ssc_ch)];
    xy_combined = [xy_combined; xy];
    
end
[r1, r2] = manual_select_population(real(log10(xy_combined)));

%% gate data all events
dead_alive_all = zeros(length(filenames),2);
for j=1:length(filenames)
    xy = [data(j).fcsdat(:,i_fsc_ch),data(j).fcsdat(:,i_ssc_ch)];
    data(j).is_dead = gate_data(real(log10(xy)), r1);
    data(j).is_alive = gate_data(real(log10(xy)), r2);
    dead_alive_all(j,1) = sum(data(j).is_dead);
    dead_alive_all(j,2) = sum(data(j).is_alive); 
end

p_dead_scatter_all = dead_alive_all(:,1)./(dead_alive_all(:,1) + dead_alive_all(:,2));




%% make fcs-ssc scatter plots
N_column = ceil(16*sqrt(length(filenames)/16/9)); % 8; % spalten
N_row = ceil(length(filenames)/N_column); % zeilen


cur_fig = figure(1); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*13 N_row*12 ], 'PaperSize', [N_column*13 N_row*12 ] );
for j=1:length(filenames)
    subplot(N_row, N_column, j)
    scatter(data(j).fcsdat(:,i_fsc_ch),data(j).fcsdat(:,i_ssc_ch), 5, (data(j).NN(:)), '.'), hold on
    

    
    dead_pgon = polyshape(10.^r1);
    plot(dead_pgon, 'FaceColor', 'none', 'EdgeColor', 'r')
    
    alive_pgon = polyshape(10.^r2);
    plot(alive_pgon, 'FaceColor', 'none', 'EdgeColor', 'g')
    
    title({sample_names{j}, [num2str(dead_alive_all(j,1)) ' dead, ' num2str(dead_alive_all(j,2)) ' alive'], [num2str(round(p_dead_scatter_all(j)*100)) '% dead'] })

        
    xlabel(data(j).fcshdr.par(i_fsc_ch).name), ylabel(data(j).fcshdr.par(i_ssc_ch).name)
    
    grid on
    caxis(NN_lim)
    set(gca,'xscale','log','yscale','log', 'XLim', scatter_lim(1:2) , 'YLim', scatter_lim(3:4) )
    
    if j == length(filenames)
        h = colorbar;
        ylabel(h, 'Nearest Neighbors')
    end
end

pause(1)
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_overview_SSC-FSC-NN.pdf']); %save figure

%%

cur_fig = figure(3); clf

subplot(2, 1, 1)
bar(p_dead_scatter_all)

set(gca, 'YLim', [0 1], 'XTick', 1:length(filenames), 'XtickLabel', {} )
grid on
ylabel('Fraction of dead cells')




subplot(2, 1, 2)
bar([dead_alive_all sum(dead_alive_all,2)]), hold on 
legend({'dead', 'alive', 'total'})
%bar(dead_alive_all(:,2)), hold on 
set(gca, 'XTick', 1:length(filenames), 'XtickLabel', sample_names )
xtickangle(30)
ylabel('Number of cells')

set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 length(filenames)*2 20 ], 'PaperSize', [length(filenames)*2 20 ] );
print(cur_fig, '-dpdf', [path_out 'fraction_dead.pdf']); %save figure



%%
disp('Done')