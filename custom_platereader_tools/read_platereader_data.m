%% Script to analyze plate reader data
close all, clear all, clc

N_channel = 3;

%% load data
[filename, pathname]=uigetfile('*.xlsx','Select data file');
path_out = pathname;
prefix_out = filename(1:end-5);

%%
tmp = xlsread([pathname filename]); % read the data only
data = tmp; %(:,3:end);

N_dp = size(data,2)/N_channel;

t = data(1, 1:N_dp)/60;
dd = data(2:end,1:N_dp);
aa = data(2:end,N_dp+1:2*N_dp);
da = data(2:end,2*N_dp+1:3*N_dp);

N_wells = size(dd,1);

%% define control wells
well_names = cell(N_wells,1);

[~, txt, ~]  = xlsread([pathname filename], 'A:A'); % read first column
[~, content, ~]  = xlsread([pathname filename], 'B:B'); % read second column

for i=1:N_wells
        well_names{i} = [num2str(i) ' ' content{i+2}];
        %well_names{i} = [txt{i+2} num2str(tmp(i))];
end
%%


N_column = 11;
N_row = ceil(N_wells/N_column);
ylim = [min([min(dd(:)) min(aa(:)) min(da(:))]) ...
    max([max(dd(:)) max(aa(:)) max(da(:)) ])];


cur_fig = figure(1); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 40 10*N_row ], 'PaperSize', [40 10*N_row] );
for i=1:N_wells
    subplot(N_row, N_column, i)
    plot(t, dd(i, :), 'g.', t, aa(i, :), 'r.', t, da(i, :), 'b.'), hold on
    hline(median(dd(i,:)), 'g');
    hline(median(aa(i,:)), 'r');
    hline(median(da(i,:)), 'b');
    set(gca, 'YLim', ylim)
    grid on
    title(well_names{i})
end

print(cur_fig, '-dpdf', [path_out filesep prefix_out '_raw.pdf']); %save figure

%%



answer = inputdlg({'Background wells:', 'Donor only:', 'Acceptor only:', 'Cells only:', 'Wells with cells'}, ...
    'Corrections',  [1 50], {'1', '2', '3', '-1', '-1'})

if str2num(answer{1})==-1
    correct_bg=0;
else
    correct_bg = 1;
    j_bg = str2num(answer{1})
end

if any(str2num(answer{2})==-1) || any(str2num(answer{3})==-1)
    correct_crosstalk = 0;
    i_donly = [];
    i_aonly = [];
else
    correct_crosstalk = 1;
    i_donly = str2num(answer{2})
    i_aonly = str2num(answer{3})
end

if any(str2num(answer{4})==-1) || any(str2num(answer{5})==-1)
    correct_cell_autofl = 0;
    i_cells_only = [];
else
    correct_cell_autofl = 1;
    i_cells_only = str2num(answer{4})
    i_cells = str2num(answer{5})
end
close(cur_fig)
%% define well names
% tmp = cell(N_wells,1);
% for i=1:N_wells
%     tmp{i} = ['Well ' num2str(i)];
% end
% 
% 
% 
% well_names = inputdlg(tmp, 'Well names',  [1 30], well_names);
%     
% titles = { 'FoB6','FoB6','D 1nM','A 1nM','D 5nM','A 5nM','D 10nM','A 10nM', ...
%     'FoB6','FoB6','LF 1nM','HF 1nM','LF 5nM','HF 5nM','LF 10nM','HF 10nM', ...
%     'FoB6','FoB6','LF 1nM','HF 1nM','LF 5nM','HF 5nM','LF 10nM','HF 10nM', ...
%     'FoB6','FoB6','D2 1nM','A2 1nM','D2 5nM','A2 5nM','D2 10nM','A2 10nM'}

%%

if correct_bg
    prefix_out = [prefix_out '_bg']
end
if correct_cell_autofl
    prefix_out = [prefix_out '_cf']
end

if correct_crosstalk
    prefix_out = [prefix_out '_ct']
end


i = 1:N_wells;
if correct_bg  
    %dd
    dd_bg = median(mean(dd(j_bg,:)));
    dd(i,:) = dd(i,:)- dd_bg;

    % aa
    aa_bg = median(mean(aa(j_bg,:)));
    aa(i,:) = aa(i,:)- aa_bg;

    % da1
    da_bg = median(mean(da(j_bg,:)));
    da(i,:) = da(i,:)- da_bg;


    

    if correct_cell_autofl
        %% determine cell fluorescence
        
        dd_bg = median(mean(dd(i_cells_only,:)));
        dd(i_cells,:) = dd(i_cells,:)- dd_bg;

        % aa
        aa_bg = median(mean(aa(i_cells_only,:)));
        aa(i_cells,:) = aa(i_cells,:)- aa_bg;

        % da1
        da_bg = mean(mean(da(i_cells_only,:)));
        da(i_cells,:) = da(i_cells,:)- da_bg;

    end
    
    
    if correct_crosstalk 
        
        % crosstalk

        leak1 = max(mean(mean(da(i_donly,:)))/mean(mean(dd(i_donly,:))),0);
        dir1 = max(mean(mean(da(i_aonly,:)))/mean(mean(aa(i_aonly,:))),0);

        da(i,:) = da(i,:)- leak1*dd(i,:)-dir1*aa(i,:); % da1
    end
end
 


%%
N_column = 8;
N_row = ceil(N_wells/N_column);
ylim = [min([min(dd(:)) min(aa(:)) min(da(:)) ]) ...
    max([max(dd(:)) max(aa(:)) max(da(:)) ])];


cur_fig = figure(1); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 40 10*N_row ], 'PaperSize', [40 10*N_row] );
for i=1:N_wells
    subplot(N_row, N_column, i)
    plot(t, dd(i, :), 'g.', t, aa(i, :), 'r.', t, da(i, :), 'b.'), hold on
    hline(mean(dd(i,:)), 'g');
    hline(mean(aa(i,:)), 'r');
    hline(mean(da(i,:)), 'b');
    set(gca, 'YLim', ylim)
    grid on
    title(well_names{i})
end

print(cur_fig, '-dpdf', [path_out filesep prefix_out '.pdf']); %save figure


%%


ylim_FRET = [-0.5 1.2];
cur_fig = figure(2); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 40 40 ], 'PaperSize', [40 40] );

i_plot = setdiff(1:N_wells, [j_bg i_donly i_aonly i_cells_only]);
for i=i_plot
    subplot(N_row, N_column, i)
    plot(t, da(i, :)./(dd(i,:)+da(i,:)), 'b.'), hold on
    hline(median(da(i,:))/(median(dd(i,:))+median(da(i,:))), 'b');
    set(gca, 'YLim', ylim_FRET)
    grid on
    title(well_names{i})
end

print(cur_fig, '-dpdf', [path_out filesep prefix_out '_FRET.pdf']); %save figure


%%


dd_mean = zeros(N_wells,3);
aa_mean = zeros(N_wells,3);
da_mean = zeros(N_wells,3);
E_mean = zeros(N_wells,3);

for i=1:N_wells
   dd_mean(i,1) = median(dd(i,:));
   dd_mean(i,2) = std(dd(i,:));
   dd_mean(i,3) = std(dd(i,:))/sqrt(length(dd(i,:)));
   
   aa_mean(i,1) = median(aa(i,:));
   aa_mean(i,2) = std(aa(i,:));
   aa_mean(i,3) = std(aa(i,:))/sqrt(length(aa(i,:)));
   
   da_mean(i,1) = median(da(i,:));
   da_mean(i,2) = std(da(i,:));
   da_mean(i,3) = std(da(i,:))/sqrt(length(da(i,:)));

   
   E_mean(i,1) = da_mean(i,1)/(dd_mean(i,1)+da_mean(i,1));
   E_mean(i,2) = sqrt(dd_mean(i,1).^2*da_mean(i,2).^2+da_mean(i,1).^2*dd_mean(i,2).^2)/(dd_mean(i,1)+da_mean(i,1)).^2;
   E_mean(i,3) = sqrt(dd_mean(i,1).^2*da_mean(i,3).^2+da_mean(i,1).^2*dd_mean(i,3).^2)/(dd_mean(i,1)+da_mean(i,1)).^2;
   
end

%
cur_fig = figure(3); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 40 10 ], 'PaperSize', [40 10] );

errorbar(1:N_wells, dd_mean(:,1), dd_mean(:,2), 'g.-'), hold on
errorbar(1:N_wells, aa_mean(:,1), aa_mean(:,2), 'r.-'), hold on
errorbar(1:N_wells, da_mean(:,1), da_mean(:,2), 'b.-'), hold on
set(gca, 'YLim', ylim, 'XLim', [0.5 N_wells+0.5], 'XTick', [1:N_wells], 'XTickLabel', well_names)
grid on
xtickangle(45)

print(cur_fig, '-dpdf', [path_out filesep prefix_out '_mean_std.pdf']); %save figure


%%
cur_fig = figure(4); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 40 10 ], 'PaperSize', [40 10] );

errorbar(i_plot, E_mean(i_plot,1), E_mean(i_plot,2), 'b.'), hold on
set(gca, 'YLim', ylim_FRET, 'XLim', [0.5 N_wells+0.5], 'XTick', [1:N_wells], 'XTickLabel', well_names)
grid on
xtickangle(45)

print(cur_fig, '-dpdf', [path_out filesep prefix_out 'FRET_mean_std.pdf']); %save figure



