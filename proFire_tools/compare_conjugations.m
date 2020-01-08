%%
close all, clear all, clc

%% select csv file and read it
n_data = 4;

data = cell(n_data,1);
for i=1:n_data
    data{i} = read_profire_csv();
    %cd(data{i}.pathname)
end
%% set manually
data_rate = 10*60; % data points/min
e = 363294; % IgG-DNA
data_rate_valve = 0.5*60; %data points / min

%fraction_times = [6.8 7.4 8.0 8.6 9.2 9.8 10.4 11 11.6 12.2 12.8 13.4 14];
%fraction = [1 2 3 4 5 6 7 8 9 10 11 12 ];

fraction_times = [7.1:0.6:14.3];
fraction = [1 2 3 4 5 6 7 8 9 10 11 12 ];

%%
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','points','PaperPosition', [0 0 3000 600], 'Position', [0 1000 1000 500]);
clf

names = cell(n_data,1);
for i=1:n_data
    plot([1:length(data{i}.Signal_numbers)]/data_rate, data{i}.Signal_numbers), hold on
    names{i} = data{i}.filename(1:end-4);
    
    
end
legend(names)

xticks = 0:3:30;

set(gca, 'XTick', xticks, 'XLim', [0 30])
grid on
%print(cur_fig, '-dpng', [data.pathname data.filename(1:end-4) '_analysis.png']); %save figure

%%

path_out = '/Users/jonasfunke/Dropbox (Personal)/Plectonic_Experiments/data_proFIRE/2019-11-08_CD33-7262-DNA-conc-screen';
prefix_out = '2019-11-08_CD33-7262-DNA-conc-screen';
cur_fig = figure(1); clf


names = cell(n_data,1);
for i=1:n_data
    plot([1:length(data{i}.Signal_numbers)]/data_rate, data{i}.Signal_numbers/abs(sum(data{i}.Signal_numbers))), hold on
    names{i} = data{i}.filename(1:end-4);
end
tmp = {'35 uM DNA', '25 uM DNA', '15 uM DNA', '5 uM DNA'};
legend(tmp, 'Interpreter', 'None')
title('CD33 with 7262, 1.1 ug/ul CD33')
xticks = 0:3:30;
xlabel('Time'), ylabel('Absorption (a.U.)')
set(gca, 'XTick', xticks, 'XLim', [2 18])
grid on
%print(cur_fig, '-dpng', [data.pathname data.filename(1:end-4) '_analysis.png']); %save figure



set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 30 20 ], 'PaperSize', [30 20] );
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_2.pdf']); %save figure



