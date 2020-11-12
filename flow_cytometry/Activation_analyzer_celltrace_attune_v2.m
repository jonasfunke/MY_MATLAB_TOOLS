%% start up
close all, clear all, clc

% set fluorescence channels
i_fsc_ch = 2; % FSC-A
i_ssc_ch = 3; % SSC-A
i_ct_ch = 4; % BL1-A for CSFE stain
i_yl1 = 6;
i_yl2 = 5;
i_rl1 = 7;
radius = 0.03;

%% load fcs data
[filenames, pathname]=uigetfile('*.fcs','Select the fcs files','MultiSelect','on');

%% create output dir
prefix_out = [ datestr(now, 'yyyy-mm-dd_HH-MM') '_activation'];
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

disp('Sample names:')
sample_names = cell(length(filenames),1);
for j=1:length(filenames)
        sample_names{j} = filenames{j}(i_found:end-4);
        disp(sample_names{j})
end

%% map back to plate

% ask if plate is used 
answer = questdlg('Did you measure one plate?', ...
	'PLate setup', ...
	'Yes','No','No');

if strcmp(answer, 'Yes')
    row = zeros(length(sample_names),1);
    column = zeros(length(sample_names),1);
    map = {'A' 1; 'B' 2; 'C' 3; 'D' 4; 'E' 5; 'F' 6; 'G' 7; 'H' 8};
    for i=1:length(sample_names)
        row(i) = find(contains(map(:,1), sample_names{i}(1)));
        column(i) = str2num(sample_names{i}(2:end));
        %disp( [sample_names{i} ' ' num2str(row(i)) ' ' num2str(column(i))] )
        data(i).row = row(i);
        data(i).column = column(i);
    end
    N_row = max(row)-min(row)+1;
    N_column = max(column)-min(column)+1;
else
    N_column = ceil(16*sqrt(length(filenames)/16/9)); % 8; % spalten
    N_row = ceil(length(filenames)/N_column); % zeilen
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


%% gate cells
ct_gate = 0.6e5; % CHANGE THIS IF NEEDED

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

%% gate cells  base on CD3 vs CD8 plots
xy_combined = zeros(0, 2);
for j=1:length(filenames)
    xy = [data(j).fcsdat(:,i_rl1),data(j).fcsdat(:,i_yl1)];
    xy_combined = [xy_combined; xy];
    
end
[r_cd8] = one_region_create_gate_2(real(log10(xy_combined)), 0.01, {['CD8 ' data(1).fcshdr.par(i_rl1).name] ['CD4 ' data(1).fcshdr.par(i_yl1).name] }, 'Select CD8+ cells');
[r_cd4] = one_region_create_gate_2(real(log10(xy_combined)), 0.01, {['CD8 ' data(1).fcshdr.par(i_rl1).name] ['CD4 ' data(1).fcshdr.par(i_yl1).name] }, 'Select CD4+ cells');


% gate cells
cd8_cd4 = zeros(length(filenames),2);
for j=1:length(filenames)
    xy = [data(j).fcsdat(:,i_rl1),data(j).fcsdat(:,i_yl1)];
    xy = real(log10(xy));
    xy(xy(:,1)<=0,1) = 1; % set inf values to one
    xy(xy(:,2)<=0,2) = 1; % set inf values to one
    data(j).cd8_positive = one_region_gate_data_2(xy, r_cd8);
    cd8_cd4(j,1) = sum(data(j).cd8_positive);
    
    xy = [data(j).fcsdat(:,i_rl1),data(j).fcsdat(:,i_yl1)];
    xy = real(log10(xy));
    xy(xy(:,1)<=0,1) = 1; % set inf values to one
    xy(xy(:,2)<=0,2) = 1; % set inf values to one
    data(j).cd4_positive = one_region_gate_data_2(xy, r_cd4);
    cd8_cd4(j,2) = sum(data(j).cd4_positive);

end


%% create scatter plots
cur_fig = figure(6); clf

for j=1:length(data)
    
    xy = [data(j).fcsdat(:,i_rl1),data(j).fcsdat(:,i_yl1)];
    xy_tmp = real([log10(xy(:,1)), log10(xy(:,2))]);
    
    
    NN = get_NN_density_fast(xy_tmp, 0.01);
    
    subplot(N_row, N_column, j)
    scatter(xy(:,1), xy(:,2), 5, NN, '.'), hold on

    poligon_tmp = polyshape(10.^r_cd4);
    plot(poligon_tmp, 'FaceColor', 'none', 'EdgeColor', 'k')
    
    poligon_tmp = polyshape(10.^r_cd8);
    plot(poligon_tmp, 'FaceColor', 'none', 'EdgeColor', 'r')
    
    title([sample_names{j} ' all'])

        
    xlabel(['CD8 ' data(j).fcshdr.par(i_rl1).name]), ylabel(['CD4 ' data(j).fcshdr.par(i_yl1).name])
    
    grid on
    if j==length(data)
        colorbar
        caxis([0 60])
        legend({'data' 'CD4+' 'CD8+'})

    end
     set(gca,'xscale','log','yscale','log', 'XLim', [1e5 1e6] , 'YLim',  [1e0 1e6])
     
     %xline(ct_gate);
     %yline(ct_gate2);
     
end

set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*14 N_row*12 ], 'PaperSize', [N_column*14 N_row*12 ] );
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_scatter_yl1-rl1.pdf']); %save figure



%% PLOT gated cells

cur_fig = figure(4); clf

for j=1:length(data)
    
    xy = [data(j).fcsdat(:,i_yl1),data(j).fcsdat(:,i_ct_ch)];
    xy_tmp = real([log10(xy(:,1)), log10(xy(:,2))]);
    
    
    NN = get_NN_density_fast(xy_tmp, 0.05);
    
    subplot(N_row, N_column, j)
    scatter(xy(:,1), xy(:,2), 5, NN, '.'), hold on

    
    title([sample_names{j} ' all'])

        
    xlabel(['CD4 ' data(j).fcshdr.par(i_yl1).name]), ylabel(['CD69/CT ' data(j).fcshdr.par(i_ct_ch).name])
    
    grid on
    if j==length(data)
        colorbar
        caxis([0 60])
    end
     set(gca,'xscale','log','yscale','log', 'XLim', [1e0 1e6] , 'YLim',  [1e0 1e6])
     
     %xline(ct_gate);
     %yline(ct_gate2);
     
end


set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*14 N_row*12 ], 'PaperSize', [N_column*14 N_row*12 ] );
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_scatter_bl1-yl1.pdf']); %save figure

%%
cur_fig = figure(5); clf

for j=1:length(data)
    
    %xy = [data(j).fcsdat(data(j).is_alive,i_rl1),data(j).fcsdat(data(j).is_alive,i_ct_ch)];
    %xy = [data(j).fcsdat(data(j).is_dead,i_rl1),data(j).fcsdat(data(j).is_dead,i_ct_ch)];
    xy = [data(j).fcsdat(:,i_rl1),data(j).fcsdat(:,i_ct_ch)];
    xy_tmp = real([log10(xy(:,1)), log10(xy(:,2))]);
    
    
    NN = get_NN_density_fast(xy_tmp, 0.01);
    
    subplot(N_row, N_column, j)
    scatter(xy(:,1), xy(:,2), 5, NN, '.'), hold on
    
   
    
    title([sample_names{j} ' all'])

        
    xlabel(['CD8 ' data(j).fcshdr.par(i_rl1).name]), ylabel(['CD69/CT ' data(j).fcshdr.par(i_ct_ch).name])
    
    grid on
    if j==length(data)
        colorbar
        caxis([0 60])
    end
     set(gca,'xscale','log','yscale','log', 'XLim', [1e5 1e6] , 'YLim',  [1e0 1e6])
     
     %xline(ct_gate);
     %yline(ct_gate2);
     
end

set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*14 N_row*12 ], 'PaperSize', [N_column*14 N_row*12 ] );
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_scatter_bl1-rl1.pdf']); %save figure


%% CD69-positive histogram
activation_gate  = 0.1e4;

activated = zeros(length(filenames),2);
cur_fig = figure(10); clf
for j=1:length(filenames)
    
    
    activated(j,1) = sum(data(j).fcsdat(data(j).cd8_positive,i_ct_ch)>activation_gate);
    activated(j,2) = sum(data(j).fcsdat(data(j).cd4_positive,i_ct_ch)>activation_gate);
    
    subplot(N_row, N_column, j)
    %histogram(real(log10(data(j).fcsdat(:,i_ct_ch))), 'DisplayStyle', 'stairs'), hold on
    histogram(real(log10(data(j).fcsdat(data(j).cd4_positive,i_ct_ch))), 'DisplayStyle', 'stairs','Normalization','pdf'), hold on
    histogram(real(log10(data(j).fcsdat(data(j).cd8_positive,i_ct_ch))), 'DisplayStyle', 'stairs', 'Normalization','pdf'), hold on
    set(gca, 'XLim', [1 6], 'YLim', [0 1])
    xline(log10(activation_gate));
    title(sample_names{j})
    if j==length(filenames)
        legend({ 'CD4+', 'CD8+'})
    end
    xlabel(['CD69, ' data(j).fcshdr.par(i_ct_ch).name])
    %ylabel('Count density')
    ylabel('PDF')
end

set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*14 N_row*12 ], 'PaperSize', [N_column*14 N_row*12 ] );
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_hist_cd69.pdf']); %save figure

%%



%% save data
close all
save([path_out prefix_out '_data.mat'])
disp('Data saved.')


%%
disp('Done')