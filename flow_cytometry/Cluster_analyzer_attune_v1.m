close all, clear all, clc

% set fluorescence channels
i_fsc_ch=2; 
i_ssc_ch=3; %FSC-A , for Cytoflex 4
i_ct_ch = 4; % BL1
i_ct_ch2 = 5; % RL1

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


%% gate cells channel 1
ct_gate = 0.8e4; % CHANGE THIS IF NEEDED

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
ct_gate2 = 0.6e4; % CHANGE THIS IF NEEDED

cur_fig = figure(1); clf
for j=1:length(filenames)
    subplot(N_row, N_column, j)
    histogram(real(log10(data(j).fcsdat(:,i_ct_ch2)))), hold on
    set(gca, 'XLim', [1 6])
    vline(log10(ct_gate2));
    data(j).is_stained2 = (data(j).fcsdat(:,i_ct_ch2)>ct_gate2);
    p_tmp = sum(data(j).is_stained2)/length(data(j).is_stained2);
    title({sample_names{j}, [ num2str(round(100*p_tmp)) '% stained cells']})
    xlabel(['CT fl, ' data(j).fcshdr.par(i_ct_ch2).name])
    ylabel('Counts')
end
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*14 N_row*12 ], 'PaperSize', [N_column*14 N_row*12 ] );
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_CT2-histogram.pdf']); %save figure


%%

cur_fig = figure(3); clf
for j=1:length(data)
    
    xy = [data(j).fcsdat(:,i_ct_ch),data(j).fcsdat(:,i_ct_ch2)];
    xy_tmp = real([log10(xy(:,1)), log10(xy(:,2))]);
    NN = get_NN_density_fast(xy_tmp, 0.05);
    
    subplot(N_row, N_column, j)
    scatter(xy(:,1), xy(:,2), 5, NN, '.'), hold on
    

    
    title(sample_names{j})

        
    xlabel(data(j).fcshdr.par(i_ct_ch).name), ylabel(data(j).fcshdr.par(i_ct_ch2).name)
    
    grid on
    colorbar
    caxis([0 60])
    set(gca,'xscale','log','yscale','log', 'XLim', [1e0 1e6] , 'YLim',  [1e0 1e6])
    
    xline(ct_gate)
    yline(ct_gate2)
end

set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 N_column*14 N_row*12 ], 'PaperSize', [N_column*14 N_row*12 ] );
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_CT1-CT2-scatter.pdf']); %save figure

%%

%%

for i=1:length(data)
    N_1(i) = sum(data(i).is_stained);
    N_2(i) = sum(data(i).is_stained2);
    N_1_2(i) = sum([data(i).is_stained & data(i).is_stained2]);
end

cur_fig = figure(4); clf

subplot(4, 1, 1:2)
bar([N_1; N_2; N_1_2]')
ylabel('CT1 positive AND CT2 positive, absolute')
set(gca, 'XTick', 1:length(sample_names), 'XTicklabel', sample_names)
legend({'CT1', 'CT2', 'CT1+CT2'})

subplot(4, 1, 3)
bar(N_1_2./N_1)
ylabel('CT1 positive AND CT2 positive / CT1')
set(gca, 'XTick', 1:length(sample_names), 'XTicklabel', sample_names, 'YLim', [0 0.5])

subplot(4, 1, 4)
bar(N_1_2./(N_1+N_2))
ylabel('CT1 positive AND CT2 positive / (CT1+CT2)')
set(gca, 'XTick', 1:length(sample_names), 'XTicklabel', sample_names, 'YLim', [0 0.1])


set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 30 50 ], 'PaperSize', [30 50] );
print(cur_fig, '-dpdf', [path_out filesep prefix_out 'counts.pdf']); %save figure


%% save data
close all
save([path_out prefix_out '_data.mat'])
disp('Data saved.')

%% export csv
file_out = [path_out prefix_out '_cluster_data.txt'];
fileID = fopen(file_out,'w');

fprintf(fileID,'Name\t');
fprintf(fileID,'N_1\t');
fprintf(fileID,'N_2\t');
fprintf(fileID,'N_1_2\t');
fprintf(fileID,'\n');

for j=1:length(filenames)
    fprintf(fileID,'%s\t', sample_names{j});
    fprintf(fileID,'%i\t',  N_1(j)); % target all
    fprintf(fileID,'%i\t', N_2(j)); % target all
    fprintf(fileID,'%i\t', N_1_2(j)); % target dead

    fprintf(fileID,'\n');
end
fclose(fileID);
disp('txt file written.')





%%
disp('Done')