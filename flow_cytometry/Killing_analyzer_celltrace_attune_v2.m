%% Killing analyzer
% v2 includes (optional) activation
%% start up
close all, clear all, clc

% set fluorescence channels
i_fsc_ch = 2; % FSC-A
i_ssc_ch = 3; % SSC-A
i_ct_ch = 4; % BL1-A for CSFE stain

i_cd4 = 6; % YL1, CD4+
i_cd8 = 7; % RL1, CD8+
i_cd69 = 5; % BL3, CD69

radius = 0.03;

%% load fcs data
[filenames, pathname]=uigetfile('*.fcs','Select the fcs files','MultiSelect','on');

%% create output dir
prefix_out = [ datestr(now, 'yyyy-mm-dd_HH-MM') '_killing-v2'];
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
answer_plate = questdlg('Did you measure one plate?', ...
	'PLate setup', ...
	'Yes','No','No');

if strcmp(answer_plate, 'Yes')
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

%% gate cells based on target-cell fluorescence
tmp = [];
for j=1:length(data)
    tmp = [tmp; data(j).fcsdat(:,i_ct_ch)];
end
ct_gate = create_gate_1d(tmp, data(j).fcshdr.par(i_ct_ch).name, 'Select gate for target cells');


cur_fig = figure(1); clf
for j=1:length(filenames)
    subplot(N_row, N_column, j)
    histogram(real(log10(data(j).fcsdat(:,i_ct_ch)))), hold on
    set(gca, 'XLim', [1 6])
    xline(log10(ct_gate));
    data(j).is_stained = (data(j).fcsdat(:,i_ct_ch)>ct_gate);
    p_tmp = sum(data(j).is_stained)/length(data(j).is_stained);
    title({sample_names{j}, [ num2str(round(100*p_tmp)) '% stained cells']})
    xlabel(['CT fl, ' data(j).fcshdr.par(i_ct_ch).name])
    ylabel('Counts')
end
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*14 N_row*12 ], 'PaperSize', [N_column*14 N_row*12 ] );
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_CT-histogram.pdf']); %save figure


%% select live and dead populations manually
xy_combined = zeros(0, 2);
for j=1:length(filenames)
    xy = [data(j).fcsdat(data(j).is_stained,i_fsc_ch),data(j).fcsdat(data(j).is_stained,i_ssc_ch)];
    xy_combined = [xy_combined; xy];
end
%[r1, r2] = manual_select_population(real(log10(xy_combined)));

r1 = create_gate_2d(xy_combined, radius, {data(1).channel_names(i_fsc_ch) data(1).channel_names(i_ssc_ch)}, 'Select dead cell population.');
r2 = create_gate_2d(xy_combined, radius, {data(1).channel_names(i_fsc_ch) data(1).channel_names(i_ssc_ch)}, 'Select alive cell population.');


%% gate data live dead target cells
dead_alive_all = zeros(length(filenames),2);
for j=1:length(filenames)
    xy = [data(j).fcsdat(data(j).is_stained,i_fsc_ch),data(j).fcsdat(data(j).is_stained,i_ssc_ch)];
    data(j).is_dead = gate_data_2d(xy, r1);
    data(j).is_alive = gate_data_2d(xy, r2);
    dead_alive_all(j,1) = sum(data(j).is_dead);
    dead_alive_all(j,2) = sum(data(j).is_alive); 
end

p_dead_scatter_all = dead_alive_all(:,1)./(dead_alive_all(:,1) + dead_alive_all(:,2));



%% calculate E:T ratio
N_target = zeros(length(filenames),1);
N_effector = zeros(length(filenames),1);
for j=1:length(filenames)
    N_effector(j) = sum(~data(j).is_stained);
    N_target(j) = sum(data(j).is_stained);
    disp(['E:T = ' num2str(round(N_effector(j)/N_target(j),1)) ', Effector: ' num2str(N_effector(j)) ', Target: ' num2str(N_target(j)) ])
end

cur_fig = figure(3); clf

subplot(2, 1, 1)
bar([N_effector N_target])
set(gca, 'XTick', 1:length(filenames), 'XtickLabel', {}, 'XLim', [0 length(filenames)+1])
grid on
ylabel('Number of counts')
legend({'Effector', 'Target'})

subplot(2, 1, 2)
bar([N_effector./N_target])

set(gca, 'XTick', 1:length(filenames), 'XtickLabel', sample_names, 'XLim', [0 length(filenames)+1], 'YLim', [0 10])
grid on
ylabel('Effector/Target')


set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 length(filenames)*2 20 ], 'PaperSize', [length(filenames)*2 20 ] );
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_effector_target_counts.pdf']); %save figure


%% make fcs-ssc scatter plots

cur_fig = figure(1); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*13 N_row*12 ], 'PaperSize', [N_column*13 N_row*12 ] );
for j=1:length(filenames)
    subplot(N_row, N_column, j)
    scatter(data(j).fcsdat(:,i_fsc_ch),data(j).fcsdat(:,i_ssc_ch), 5, (data(j).NN(:)), '.'), hold on
    

    
    dead_pgon = polyshape(r1);
    plot(dead_pgon, 'FaceColor', 'none', 'EdgeColor', 'r')
    
    alive_pgon = polyshape(r2);
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

%
cur_fig = figure(1); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*13 N_row*12 ], 'PaperSize', [N_column*13 N_row*12 ] );
for j=1:length(filenames)
    
    xy = [data(j).fcsdat(data(j).is_stained,i_fsc_ch),data(j).fcsdat(data(j).is_stained,i_ssc_ch)];
    xy_tmp = real([log10(xy(:,1)), log10(xy(:,2))]);
    NN = get_NN_density_fast(xy_tmp, radius);
    
    subplot(N_row, N_column, j)
    scatter(xy(:,1), xy(:,2), 5, NN, '.'), hold on
    

    
    dead_pgon = polyshape(r1);
    plot(dead_pgon, 'FaceColor', 'none', 'EdgeColor', 'r')
    
    alive_pgon = polyshape(r2);
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

% plot pbmc only
cur_fig = figure(1); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*13 N_row*12 ], 'PaperSize', [N_column*13 N_row*12 ] );
for j=1:length(filenames)
    
    xy = [data(j).fcsdat(~data(j).is_stained,i_fsc_ch),data(j).fcsdat(~data(j).is_stained,i_ssc_ch)];
    xy_tmp = real([log10(xy(:,1)), log10(xy(:,2))]);
    NN = get_NN_density_fast(xy_tmp, radius);
    
    subplot(N_row, N_column, j)
    scatter(xy(:,1), xy(:,2), 5, NN, '.'), hold on
    

    
    dead_pgon = polyshape(r1);
    plot(dead_pgon, 'FaceColor', 'none', 'EdgeColor', 'r')
    
    alive_pgon = polyshape(r2);
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
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_overview_SSC-FSC-NN_not-stained.pdf']); %save figure

%%

cur_fig = figure(3); clf

subplot(2, 1, 1)
bar(p_dead_scatter_all)

set(gca, 'YLim', [0 1], 'XTick', 1:length(filenames), 'XtickLabel', {}, 'XLim', [0 length(filenames)+1])
grid on
ylabel('Fraction of dead cells')




subplot(2, 1, 2)
bar([dead_alive_all sum(dead_alive_all,2)]), hold on 
legend({'dead', 'alive', 'total'})
%bar(dead_alive_all(:,2)), hold on 
set(gca, 'XTick', 1:length(filenames), 'XtickLabel', sample_names, 'XLim', [0 length(filenames)+1] )
xtickangle(30)
ylabel('Number of cells')
grid on

set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 length(filenames)*2 20 ], 'PaperSize', [length(filenames)*2 20 ] );
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_fraction_dead.pdf']); %save figure


%%

if answer_plate
    dead_alive_plate = zeros(8, 12, 2);

    for i=1:length(filenames)
        dead_alive_plate(row(i), column(i), :) = dead_alive_all(i,:);    
    end


    cur_fig = figure(4); clf

    subplot(4, 3, [1:9])
    imagesc(dead_alive_plate(:,:,1)./(dead_alive_plate(:,:,1) + dead_alive_plate(:,:,2)), [0 1]), axis image, colorbar
    set(gca, 'Xtick', 1:12, 'YTick', [1:8], 'Yticklabel', {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'})
    title('Fraction of dead target cells')

    subplot(4, 3, 10)
    imagesc(dead_alive_plate(:,:,2)), axis image, colorbar
    set(gca, 'Xtick', 1:12, 'YTick', [1:8], 'Yticklabel', {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'})
    title('Alive target cells')


    subplot(4, 3, 11)
    imagesc(dead_alive_plate(:,:,1)), axis image, colorbar
    set(gca, 'Xtick', 1:12, 'YTick', [1:8], 'Yticklabel', {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'})
    title('Dead target cells')
    
     subplot(4, 3, 12)
    imagesc(dead_alive_plate(:,:,1)+dead_alive_plate(:,:,2)), axis image, colorbar
    set(gca, 'Xtick', 1:12, 'YTick', [1:8], 'Yticklabel', {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'})
    title('Total target cells')

    set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
        'PaperPosition', [0 0 30 30 ], ...
        'PaperSize', [30 30] );

    print(cur_fig, '-dpdf', [path_out filesep prefix_out '_killing_image.pdf']); %save figure
end

%% all events

N_sample = length(filenames);

N_counts = zeros(N_sample,1);

for j=1:length(filenames)
    N_counts(j,1) = length(data(j).fcsdat(:,1));
end





%% calculate activation
close all
answer_activation = questdlg('Analyze T-cell activation?', ...
	'Activation', ...
	'Yes','No','No');
bool_activation = strcmp(answer_activation, 'Yes');

if bool_activation

    % gate cells base on CD3 vs CD8 plots
    xy_combined = zeros(0, 2);
    for j=1:length(filenames)
        xy = [data(j).fcsdat(:,i_cd8),data(j).fcsdat(:,i_cd4)];
        xy_combined = [xy_combined; xy];
    end
    
    %[r_cd8] = one_region_create_gate_2(real(log10(xy_combined)), 0.01, {['CD8 ' data(1).fcshdr.par(i_rl1).name] ['CD4 ' data(1).fcshdr.par(i_yl1).name] }, 'Select CD8+ cells');
    %[r_cd4] = one_region_create_gate_2(real(log10(xy_combined)), 0.01, {['CD8 ' data(1).fcshdr.par(i_rl1).name] ['CD4 ' data(1).fcshdr.par(i_yl1).name] }, 'Select CD4+ cells');

    cd4_gate = create_gate_1d(xy_combined(:,2), ['CD4 ' data(1).fcshdr.par(i_cd4).name], 'Select gate for CD4+ cells');
    cd8_gate = create_gate_1d(xy_combined(:,1), ['CD8 ' data(1).fcshdr.par(i_cd8).name], 'Select gate for CD8+ cells');

    % gate cells
    cd8_cd4 = zeros(length(filenames),2);
    for j=1:length(filenames)
        %xy = [data(j).fcsdat(:,i_rl1),data(j).fcsdat(:,i_yl1)];
        %xy = real(log10(xy));
        %xy(xy(:,1)<=0,1) = 1; % set inf values to one
        %xy(xy(:,2)<=0,2) = 1; % set inf values to one
        %data(j).cd8_positive = one_region_gate_data_2(xy, r_cd8);
        data(j).cd8_positive = data(j).fcsdat(:,i_cd8)>cd8_gate;
        cd8_cd4(j,1) = sum(data(j).cd8_positive);

        %xy = [data(j).fcsdat(:,i_rl1),data(j).fcsdat(:,i_yl1)];
        %xy = real(log10(xy));
        %xy(xy(:,1)<=0,1) = 1; % set inf values to one
        %xy(xy(:,2)<=0,2) = 1; % set inf values to one
        %data(j).cd4_positive = one_region_gate_data_2(xy, r_cd4);
        data(j).cd4_positive = data(j).fcsdat(:,i_cd4)>cd4_gate;
        cd8_cd4(j,2) = sum(data(j).cd4_positive);

    end



    %% create scatter plots for activation
    cur_fig = figure(6); clf

    for j=1:length(data)

        xy = [data(j).fcsdat(:,i_cd8),data(j).fcsdat(:,i_cd4)];
        xy_tmp = real([log10(xy(:,1)), log10(xy(:,2))]);


        NN = get_NN_density_fast(xy_tmp, 0.05);

        subplot(N_row, N_column, j)
        scatter(xy(:,1), xy(:,2), 5, NN, '.'), hold on

        %poligon_tmp = polyshape(10.^r_cd4);
        %plot(poligon_tmp, 'FaceColor', 'none', 'EdgeColor', 'k')
        yline(cd4_gate, 'r');
        %poligon_tmp = polyshape(10.^r_cd8);
        %plot(poligon_tmp, 'FaceColor', 'none', 'EdgeColor', 'r')
        xline(cd8_gate, 'r');
        title([sample_names{j} ' all'])


        xlabel(['CD8 ' data(j).fcshdr.par(i_cd8).name]), ylabel(['CD4 ' data(j).fcshdr.par(i_cd4).name])

        grid on

        caxis([0 60])
        if j==length(data)
            colorbar
            legend({'data' 'CD4+' 'CD8+'})

        end
         set(gca,'xscale','log','yscale','log', 'XLim', [1e0 1e6] , 'YLim',  [1e0 1e6])

         %xline(ct_gate);
         %yline(ct_gate2);

    end

    set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
        'PaperPosition', [0 0 N_column*14 N_row*12 ], 'PaperSize', [N_column*14 N_row*12 ] );
    print(cur_fig, '-dpdf', [path_out filesep prefix_out '_scatter_CD4-CD4.pdf']); %save figure



    %% Plot gated cells

    cur_fig = figure(4); clf

    for j=1:length(data)

        xy = [data(j).fcsdat(:,i_cd4),data(j).fcsdat(:,i_cd69)];
        xy_tmp = real([log10(xy(:,1)), log10(xy(:,2))]);


        NN = get_NN_density_fast(xy_tmp, 0.05);

        subplot(N_row, N_column, j)
        scatter(xy(:,1), xy(:,2), 5, NN, '.'), hold on
        xline(cd4_gate, 'r');


        title([sample_names{j} ' all'])


        xlabel(['CD4 ' data(j).fcshdr.par(i_cd4).name]), ylabel(['CD69 ' data(j).fcshdr.par(i_cd69).name])
        caxis([0 60])

        grid on
        if j==length(data)
            colorbar
    %        caxis([0 60])
        end

         set(gca,'xscale','log','yscale','log', 'XLim', [1e0 1e6] , 'YLim',  [1e0 1e6])

         %xline(ct_gate);
         %yline(ct_gate2);

    end


    set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
        'PaperPosition', [0 0 N_column*14 N_row*12 ], 'PaperSize', [N_column*14 N_row*12 ] );
    print(cur_fig, '-dpdf', [path_out filesep prefix_out '_scatter_CD69-CD4.pdf']); %save figure

    %%
    cur_fig = figure(5); clf

    for j=1:length(data)

        %xy = [data(j).fcsdat(data(j).is_alive,i_rl1),data(j).fcsdat(data(j).is_alive,i_ct_ch)];
        %xy = [data(j).fcsdat(data(j).is_dead,i_rl1),data(j).fcsdat(data(j).is_dead,i_ct_ch)];
        xy = [data(j).fcsdat(:,i_cd8),data(j).fcsdat(:,i_cd69)];
        xy_tmp = real([log10(xy(:,1)), log10(xy(:,2))]);


        NN = get_NN_density_fast(xy_tmp, 0.05);

        subplot(N_row, N_column, j)
        scatter(xy(:,1), xy(:,2), 5, NN, '.'), hold on
        xline(cd8_gate, 'r');

        % plot cd8+ cells
        %xy = [data(j).fcsdat(data(j).cd8_positive,i_rl1),data(j).fcsdat(data(j).cd8_positive,i_ct_ch)];
        %scatter(xy(:,1), xy(:,2), 5, 'r.'), hold on




        title([sample_names{j} ' all'])


        xlabel(['CD8 ' data(j).fcshdr.par(i_cd8).name]), ylabel(['CD69 ' data(j).fcshdr.par(i_cd69).name])
        caxis([0 60])

        grid on
        if j==length(data)
            colorbar
    %        caxis([0 40])
        end

         set(gca,'xscale','log','yscale','log', 'XLim', [1e0 1e6] , 'YLim',  [1e0 1e6])

         %xline(ct_gate);
         %yline(ct_gate2);

    end

    set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
        'PaperPosition', [0 0 N_column*14 N_row*12 ], 'PaperSize', [N_column*14 N_row*12 ] );
    print(cur_fig, '-dpdf', [path_out filesep prefix_out '_scatter_CD69-CD8.pdf']); %save figure


    %% CD69-positive histogram
    %activation_gate_cd4  = 0.5e4;
    %activation_gate_cd8  = 0.5e4;

    tmp1 = [];
    tmp2 = [];
    for j=1:length(filenames)
        tmp1 = [tmp1; data(j).fcsdat(data(j).cd4_positive,i_cd69)];
        tmp2 = [tmp2; data(j).fcsdat(data(j).cd8_positive,i_cd69)];
    end
    activation_gate_cd4 = create_gate_1d(tmp1, data(j).fcshdr.par(i_cd69).name, 'Select gate for CD4+ cells');
    activation_gate_cd8 = create_gate_1d(tmp2, data(j).fcshdr.par(i_cd69).name, 'Select gate for CD8+ cells');



    cc = lines(2);

    activated = zeros(length(filenames),2);
    cur_fig = figure(10); clf
    for j=1:length(filenames)


        activated(j,2) = sum(data(j).fcsdat(data(j).cd8_positive,i_cd69)>activation_gate_cd8);
        activated(j,1) = sum(data(j).fcsdat(data(j).cd4_positive,i_cd69)>activation_gate_cd4);

        subplot(N_row, N_column, j)
        %histogram(real(log10(data(j).fcsdat(:,i_ct_ch))), 'DisplayStyle', 'stairs'), hold on
        histogram(real(log10(data(j).fcsdat(data(j).cd4_positive,i_cd69))), 'DisplayStyle', 'stairs', 'Normalization','pdf'), hold on
        histogram(real(log10(data(j).fcsdat(data(j).cd8_positive,i_cd69))), 'DisplayStyle', 'stairs', 'Normalization','pdf'), hold on
        set(gca, 'XLim', [1 6], 'YLim', [0 1])
        xline(log10(activation_gate_cd4), 'Color', cc(1,:));
        xline(log10(activation_gate_cd8), 'Color', cc(2,:));
        
        title({sample_names{j}, ...
            [num2str(round(100*activated(j,1)/sum(data(j).cd4_positive))) '% activated CD4'], ...
            [num2str(round(100*activated(j,2)/sum(data(j).cd8_positive))) '% activated CD8'] })
        if j==length(filenames)
            legend({ 'CD4+', 'CD8+'})
        end
        xlabel(['CD69, ' data(j).fcshdr.par(i_cd69).name])
        %ylabel('Count density')
        ylabel('PDF')
    end

    set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
        'PaperPosition', [0 0 N_column*14 N_row*12 ], 'PaperSize', [N_column*14 N_row*12 ] );
    print(cur_fig, '-dpdf', [path_out filesep prefix_out '_hist_cd69.pdf']); %save figure

    %%

    if answer_plate
        plate_activated_cd4 = zeros(8, 12, 2);
        plate_activated_cd8 = zeros(8, 12, 2);
        for i=1:length(filenames)
            plate_activated_cd4(row(i), column(i), 1) = activated(i,1);   % activated cd4 cells
            plate_activated_cd4(row(i), column(i), 2) = sum(data(i).cd4_positive); % number of cd4 cells

            plate_activated_cd8(row(i), column(i), 1) = activated(i,2);   % activated cd4 cells
            plate_activated_cd8(row(i), column(i), 2) = sum(data(i).cd8_positive); % number of cd4 cells

        end


        cur_fig = figure(4); clf

        subplot(2, 2, 1)
        imagesc(plate_activated_cd4(:,:,1)), axis image, colorbar
        set(gca, 'Xtick', 1:12, 'YTick', [1:8], 'Yticklabel', {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'})
        title(['Number activated CD4+ ' data(j).fcshdr.par(i_cd4).name ' cells'])

        subplot(2, 2, 2)
        imagesc(plate_activated_cd4(:,:,1)./plate_activated_cd4(:,:,2), [0 1]), axis image, colorbar
        set(gca, 'Xtick', 1:12, 'YTick', [1:8], 'Yticklabel', {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'})
        title(['Fraction activated CD4+ ' data(j).fcshdr.par(i_cd4).name ' cells'])



        subplot(2, 2, 3)
        imagesc(plate_activated_cd8(:,:,1)), axis image, colorbar
        set(gca, 'Xtick', 1:12, 'YTick', [1:8], 'Yticklabel', {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'})
        title(['Number activated CD8+ ' data(j).fcshdr.par(i_cd8).name ' cells'])

        subplot(2, 2, 4)
        imagesc(plate_activated_cd8(:,:,1)./plate_activated_cd8(:,:,2), [0 1]), axis image, colorbar
        set(gca, 'Xtick', 1:12, 'YTick', [1:8], 'Yticklabel', {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'})
        title(['Fraction activated CD8+ ' data(j).fcshdr.par(i_cd8).name ' cells'])


        set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
            'PaperPosition', [0 0 40 20 ], ...
            'PaperSize', [40 20] );

        print(cur_fig, '-dpdf', [path_out filesep prefix_out '_activation_image.pdf']); %save figure
    end

end
%% calculate fraction of cd8+ cells
% tmp = zeros(8,12);
% tmp2 = zeros(8,12);
% 
% for i=1:length(filenames)
% 
%     tmp(row(i), column(i)) = sum(data(i).cd8_positive)./length(data(i).cd8_positive);
%     tmp2(row(i), column(i)) = sum(data(i).is_stained)./length(data(i).is_stained);
% 
% end
%% save data
file_out = [path_out prefix_out '_killing_data.txt'];
fileID = fopen(file_out,'w');

fprintf(fileID,'Name\t');
fprintf(fileID,'N_all\t');
fprintf(fileID,'N_target\t');
fprintf(fileID,'N_target_dead\t');
fprintf(fileID,'N_target_alive\t');
if bool_activation
    fprintf(fileID,'N_cd4_total\t');
    fprintf(fileID,'N_cd4_acivated\t');
    fprintf(fileID,'N_cd8_total\t');
    fprintf(fileID,'N_cd8_acivated\t');
end

fprintf(fileID,'\n');

for j=1:length(filenames)
    fprintf(fileID,'%s\t', sample_names{j});
    fprintf(fileID,'%i\t',  N_counts(j,1)); % target all
    fprintf(fileID,'%i\t', sum(dead_alive_all(j,:))); % target all
    fprintf(fileID,'%i\t', dead_alive_all(j,1)); % target dead
    fprintf(fileID,'%i\t', dead_alive_all(j,2)); % target alive
    if bool_activation
        fprintf(fileID,'%i\t',sum(data(j).cd4_positive)); % cd4 total
        fprintf(fileID,'%i\t', activated(j,1)); % cd4 activated
        fprintf(fileID,'%i\t',sum(data(j).cd8_positive)); % cd8 total
        fprintf(fileID,'%i\t', activated(j,2)); % cd8 activated
    end
    fprintf(fileID,'\n');
end
fclose(fileID);
disp('txt file written.')



%%
close all
clear('tmp', 'tmp1', 'tmp2', 'xy', 'xy_combined', 'xy_tmp', 'NN')
save([path_out prefix_out '_data.mat'])
disp('Data saved.')






%%
disp('Done')