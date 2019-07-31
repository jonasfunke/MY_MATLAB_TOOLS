%%
close all, clear all, clc

%% select csv file and read it

data = read_profire_csv();

%% set manually

data_rate = 10*60; % data points/min
e = 363294; % IgG-DNA
data_rate_valve = 0.5*60; %data points / min

% custom_BL_IgG_collect
%fraction_times = [3.0:0.6:4.8 8.6:0.6:14];
%fraction = [1 2 3 13 4 5 6 7 8 9 10 11 12 ];

% custom_49_IgG-collect
%fraction_times = [2.55 3.2 3.8 4.4 9.72:0.66:15.66];
%fraction = [1 2 3 13 4 5 6 7 8 9 10 11 12 ];

% custom-BK-IgG-26
fraction_times = [6.8 7.4 8.0 8.6 9.2 9.8 10.4 11 11.6 12.2 12.8 13.4 14];
fraction = [1 2 3 4 5 6 7 8 9 10 11 12 ];


%%

cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','points','PaperPosition', [0 0 3000 600], 'Position', [0 1000 1000 500]);
clf
plot([1:length(data.Signal_numbers)]/data_rate, data.Signal_numbers), hold on

for j=1:length(fraction_times)-1
    vline(fraction_times(j), 'k:');
    
    index(1) = round(data_rate*fraction_times(j));
    index(2) = round(data_rate*fraction_times(j+1));
    
    
    avg_abs = mean(data.Signal_numbers(index(1):index(2)));

    dt = (index(2)-index(1))/data_rate; % time for fraction
    dV = dt*1*1e3; % Volume in fraciton (in ul) given 1ml/min flow rate

    concentration = 1e9*avg_abs/e/dV;
    amount = 1e6*avg_abs/e; %pmol

    text(fraction_times(j)+0.1*(fraction_times(j+1)-fraction_times(j)), 0.05, { ['F' num2str(fraction(j))],[num2str(round(concentration)) ' nM'], [num2str(round(amount)) ' pmol']})

end
vline(fraction_times(length(fraction_times)), 'k:');
title(['Using 363294/M/cm extinction coefficient'])
set(gca, 'XTick', 0:3:27, 'XLim', [0 30])%, 'YLim', [-10 130])
grid on
print(cur_fig, '-dpng', [data.pathname data.filename(1:end-4) '_analysis.png']); %save figure


