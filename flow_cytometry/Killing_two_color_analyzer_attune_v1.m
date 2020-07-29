close all, clear all, clc

% set fluorescence channels
i_fsc_ch=2; 
i_ssc_ch=3; %FSC-A , for Cytoflex 4
i_ct_ch = 4; % BL1
i_ct2_ch = 5; % RL1

radius = 0.03;

%%
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
    data(i) = load_fcs_data_attune_v2(pathname, filenames{i}, radius);
end

%% create sample names
i=1;
while i<length(filenames{1})
    pattern = filenames{1}(1:end-i);
    if all(startsWith(filenames, pattern))
        i_found=length(filenames{1})-i+1;
        i=length(filenames{1}); % stop
    end
    i = i+1;
end

sample_names = cell(length(filenames),1);
for j=1:length(filenames)
        sample_names{j} = filenames{j}(i_found:end-4);
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

N_column = ceil(16*sqrt(length(filenames)/16/9)); % 8; % spalten
N_row = ceil(length(filenames)/N_column); % zeilen

%% gate cells
ct_gate = 1e4; % CHANGE THIS IF NEEDED

cur_fig = figure(1); clf
for j=1:length(filenames)
    subplot(N_row, N_column, j)
    histogram(real(log10(data(j).fcsdat(:,i_ct_ch)))), hold on
    set(gca, 'XLim', [1 6])
    vline(log10(ct_gate));
    data(j).is_stained = (data(j).fcsdat(:,i_ct_ch)>ct_gate);
    p_tmp = sum(data(j).is_stained)/length(data(j).is_stained);
    title({sample_names{j}, [ num2str(round(100*p_tmp)) '% stained cells']})
    xlabel(['CT fl, ' data(j).fcshdr.par(i_ct_ch).name])
    ylabel('Counts')
end
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*14 N_row*12 ], 'PaperSize', [N_column*14 N_row*12 ] );
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_CT-histogram.pdf']); %save figure


%% gate cells 2
ct_gate_2 = 0.7e4; % CHANGE THIS IF NEEDED

cur_fig = figure(1); clf
for j=1:length(filenames)
    subplot(N_row, N_column, j)
    histogram(real(log10(data(j).fcsdat(:,i_ct2_ch)))), hold on
    set(gca, 'XLim', [1 6])
    vline(log10(ct_gate_2));
    data(j).is_stained2 = (data(j).fcsdat(:,i_ct2_ch)>ct_gate_2);
    p_tmp = sum(data(j).is_stained2)/length(data(j).is_stained2);
    title({sample_names{j}, [ num2str(round(100*p_tmp)) '% stained cells']})
    xlabel(['CT fl, ' data(j).fcshdr.par(i_ct2_ch).name])
    ylabel('Counts')
end
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*14 N_row*12 ], 'PaperSize', [N_column*14 N_row*12 ] );
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_CT2-histogram.pdf']); %save figure


%%
xy_combined = zeros(0, 2);
for j=1:length(filenames)
    xy = [data(j).fcsdat(data(j).is_stained,i_fsc_ch),data(j).fcsdat(data(j).is_stained,i_ssc_ch)];
    xy_combined = [xy_combined; xy];
    
end
[r1, r2] = manual_select_population(real(log10(xy_combined)));

% gate data all events
dead_alive_all = zeros(length(filenames),2);
for j=1:length(filenames)
    xy = [data(j).fcsdat(data(j).is_stained,i_fsc_ch),data(j).fcsdat(data(j).is_stained,i_ssc_ch)];
    data(j).is_dead = gate_data(real(log10(xy)), r1);
    data(j).is_alive = gate_data(real(log10(xy)), r2);
    dead_alive_all(j,1) = sum(data(j).is_dead);
    dead_alive_all(j,2) = sum(data(j).is_alive); 
end

p_dead_scatter_all = dead_alive_all(:,1)./(dead_alive_all(:,1) + dead_alive_all(:,2));


%% 2nd channel
xy_combined = zeros(0, 2);
for j=1:length(filenames)
    xy = [data(j).fcsdat(data(j).is_stained2,i_fsc_ch),data(j).fcsdat(data(j).is_stained2,i_ssc_ch)];
    xy_combined = [xy_combined; xy];
    
end
[r1_2, r2_2] = manual_select_population(real(log10(xy_combined)));

%% gate data all events
dead_alive_all2 = zeros(length(filenames),2);
for j=1:length(filenames)
    xy = [data(j).fcsdat(data(j).is_stained2,i_fsc_ch),data(j).fcsdat(data(j).is_stained2,i_ssc_ch)];
    data(j).is_dead2 = gate_data(real(log10(xy)), r1_2);
    data(j).is_alive2 = gate_data(real(log10(xy)), r2_2);
    dead_alive_all2(j,1) = sum(data(j).is_dead2);
    dead_alive_all2(j,2) = sum(data(j).is_alive2); 
end

p_dead_scatter_all2 = dead_alive_all2(:,1)./(dead_alive_all2(:,1) + dead_alive_all2(:,2));




%% make fcs-ssc scatter plots



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
cur_fig = figure(1); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*13 N_row*12 ], 'PaperSize', [N_column*13 N_row*12 ] );
for j=1:length(filenames)
    
    xy = [data(j).fcsdat(data(j).is_stained,i_fsc_ch),data(j).fcsdat(data(j).is_stained,i_ssc_ch)];
    xy_tmp = real([log10(xy(:,1)), log10(xy(:,2))]);
    NN = get_NN_density_fast(xy_tmp, radius);
    
    subplot(N_row, N_column, j)
    scatter(xy(:,1), xy(:,2), 5, NN, '.'), hold on
    

    
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
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_overview_SSC-FSC-NN_stained.pdf']); %save figure
%% 2
cur_fig = figure(1); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*13 N_row*12 ], 'PaperSize', [N_column*13 N_row*12 ] );
for j=1:length(filenames)
    
    xy = [data(j).fcsdat(data(j).is_stained2,i_fsc_ch),data(j).fcsdat(data(j).is_stained2,i_ssc_ch)];
    xy_tmp = real([log10(xy(:,1)), log10(xy(:,2))]);
    NN = get_NN_density_fast(xy_tmp, radius);
    
    subplot(N_row, N_column, j)
    scatter(xy(:,1), xy(:,2), 5, NN, '.'), hold on
    

    
    dead_pgon = polyshape(10.^r1_2);
    plot(dead_pgon, 'FaceColor', 'none', 'EdgeColor', 'r')
    
    alive_pgon = polyshape(10.^r2_2);
    plot(alive_pgon, 'FaceColor', 'none', 'EdgeColor', 'g')
    
    title({sample_names{j}, [num2str(dead_alive_all2(j,1)) ' dead, ' num2str(dead_alive_all2(j,2)) ' alive'], [num2str(round(p_dead_scatter_all2(j)*100)) '% dead'] })

        
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
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_overview_SSC-FSC-NN_stained2.pdf']); %save figure



%% plot pbmc only
cur_fig = figure(1); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*13 N_row*12 ], 'PaperSize', [N_column*13 N_row*12 ] );
for j=1:length(filenames)
    i_use = ~data(j).is_stained & ~data(j).is_stained2; 
    xy = [data(j).fcsdat(i_use,i_fsc_ch),data(j).fcsdat(i_use,i_ssc_ch)];
    xy_tmp = real([log10(xy(:,1)), log10(xy(:,2))]);
    NN = get_NN_density_fast(xy_tmp, radius);
    
    subplot(N_row, N_column, j)
    scatter(xy(:,1), xy(:,2), 5, NN, '.'), hold on
    
    
    title({sample_names{j} })

        
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
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_overview_SSC-FSC-NN_not-stained.pdf']); %save figure

%%

cur_fig = figure(3); clf

subplot(3, 1, 1)
bar([p_dead_scatter_all p_dead_scatter_all2])
legend({'BL1', 'RL1'})
set(gca, 'YLim', [0 1], 'XTick', 1:length(filenames), 'XtickLabel', {}, 'XLim', [0 length(filenames)+1])
grid on
ylabel('Fraction of dead cells')




subplot(3, 1, 2)
bar([dead_alive_all sum(dead_alive_all,2)]), hold on 
legend({'dead', 'alive', 'total'})
%bar(dead_alive_all(:,2)), hold on 
set(gca, 'XTick', 1:length(filenames), 'XtickLabel', sample_names, 'XLim', [0 length(filenames)+1] )
xtickangle(30)
ylabel('Number of cells, BL1')
grid on


subplot(3, 1, 3)
bar([dead_alive_all2 sum(dead_alive_all2,2)]), hold on 
legend({'dead', 'alive', 'total'})
%bar(dead_alive_all(:,2)), hold on 
set(gca, 'XTick', 1:length(filenames), 'XtickLabel', sample_names, 'XLim', [0 length(filenames)+1] )
xtickangle(30)
ylabel('Number of cells, RL1')
grid on


set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 length(filenames)*2 20 ], 'PaperSize', [length(filenames)*2 20 ] );
print(cur_fig, '-dpdf', [path_out 'fraction_dead.pdf']); %save figure



%% all events

N_sample = length(filenames);

N_counts = zeros(N_sample,1);

for j=1:length(filenames)
    N_counts(j,1) = length(data(j).fcsdat(:,1));
end

%% save data
close all
save([path_out prefix_out '_data.mat'])
disp('Data saved.')

%% export csv
file_out = [path_out prefix_out '_killing_data.txt'];
fileID = fopen(file_out,'w');

fprintf(fileID,'Name\t');
fprintf(fileID,'N_all\t');
fprintf(fileID,'N_target\t');
fprintf(fileID,'N_target_dead\t');
fprintf(fileID,'N_target_alive\t');
fprintf(fileID,'\n');

for j=1:length(filenames)
    fprintf(fileID,'%s\t', sample_names{j});
    fprintf(fileID,'%i\t',  N_counts(j,1)); % target all
    fprintf(fileID,'%i\t', sum(dead_alive_all(j,:))); % target all
    fprintf(fileID,'%i\t', dead_alive_all(j,1)); % target dead
    fprintf(fileID,'%i\t', dead_alive_all(j,2)); % target alive

    fprintf(fileID,'\n');
end
fclose(fileID);
disp('txt file written.')





%%
disp('Done')