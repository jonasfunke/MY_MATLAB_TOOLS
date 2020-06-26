%%
close all, clear all, clc

%% select csv file and read it
data = read_profire_csv();

%% set manually

data_rate = 10*60; % data points/min
%e = 363294; % IgG-DNA
%e = 136700 + 11710; % IL2 short modifier + extinction at 280
e = 225800+ 11710; % IL2 long modifier + extinction at 280

data_rate_valve = 0.5*60; %data points / min
%%

% e_IgG = 141894;
% e_7020 = 221400*15/26;
% e = e_IgG + e_7020;

% custom-BK-IgG-26_collect
% fraction_times = [3.0:0.6:4.8 8.6:0.6:14];
% fraction = [1 2 3 13 4 5 6 7 8 9 10 11 12 ];

% custom-IL2-short
% fraction_times = [11.4:0.6:18.6];
% fraction = [1:12];

%custom-IL2-longv2
fraction_times = [8.4:0.6:15.6];
fraction = [1:12];

% custom_49_IgG-collect
%fraction_times = [2.55 3.2 3.8 4.4 9.72:0.66:15.66];
%fraction = [1 2 3 13 4 5 6 7 8 9 10 11 12 ];

% %custom-BK-IgG-26
% fraction_times = [6.8 7.4 8.0 8.6 9.2 9.8 10.4 11 11.6 12.2 12.8 13.4 14];
% fraction = [1 2 3 4 5 6 7 8 9 10 11 12 ];

% % custom-BK-IgG-26-all-peaks
% fraction_times = [9.8 10.4 11 11.6 12.2 12.8 13.4 14:0.6:17];
% fraction = [1 2 3 4 5 6 7 8 9 10 11 12 ];

%proFIRE19_32-41bases
% fraction_times = [7.2:0.65:15];
% fraction = [1 2 3 4 5 6 7 8 9 10 11 12 ];

%proFIRE16_29-31bases
% fraction_times = [6.0:0.6:13.2];
% fraction = [1 2 3 4 5 6 7 8 9 10 11 12 ];

%proFIRE10_15-19bases
% fraction_times = [6.0:0.6:13.2];
% fraction = [1 2 3 4 5 6 7 8 9 10 11 12 ];

% custom-BK-IgG-29-31
% fraction_times = [7.1:0.6:14.3];
% fraction = [1 2 3 4 5 6 7 8 9 10 11 12 ];


%% determine fraction_times and fractions
%  figure(1), clf
% tmp = 1:length(data.Position_numbers);
% % 
% % plot( data.Position_numbers)
% plot(tmp/data_rate_valve, data.Position_numbers)
% %%
% fraction = unique(data.Position_numbers);
% fraction = fraction(1:end-1); %remove fraction 13
% fraction_times = zeros(1, length(fraction)+1);
% for j=fraction
%     index = tmp(data.Position_numbers==j); % index of fractions
%     selection = index-index(1)-cumsum(ones(1,length(index)))+1; % select only thefirst fraction
%     index = index(selection==0);
%     fraction_times(j) = (index(1)-1)/data_rate_valve;
% end
% fraction_times(end) = (index(end))/data_rate_valve
% %fraction
%%

cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 80 12 ], 'PaperSize', [80 12] );
clf
plot([1:length(data.Signal_numbers)]/data_rate, data.Signal_numbers), hold on
xlim = [0 max([1:length(data.Signal_numbers)]/data_rate)];
ylim = [min(data.Signal_numbers) 1.1*max(data.Signal_numbers)];
xtick = 0:2:ylim(2);
for j=1:length(fraction_times)-1
    vline(fraction_times(j), 'k-');
    
    index(1) = round(data_rate*fraction_times(j));
    index(2) = round(data_rate*fraction_times(j+1));
    
    
    avg_abs = mean(data.Signal_numbers(index(1):index(2)));

    dt = (index(2)-index(1))/data_rate; % time for fraction
    dV = dt*1*1e3; % Volume in fraciton (in ul) given 1ml/min flow rate

    concentration = 1e9*avg_abs/e/dV;
    amount = 1e6*avg_abs/e; %pmol

    text(fraction_times(j)+0.1*(fraction_times(j+1)-fraction_times(j)), (ylim(2)-ylim(1))/2, { ['F' num2str(fraction(j))],[num2str(round(concentration)) ' nM'], [num2str(round(amount)) ' pmol']})

end
vline(fraction_times(length(fraction_times)), 'k-');
title([ data.filename(1:end-4) ', extinction coefficient=' num2str(e) ' /M/cm'])
set(gca, 'XLim', xlim, 'YLim', ylim, 'XTick', xtick)
grid on

%print(cur_fig, '-dpng', [data.pathname data.filename(1:end-4) '_analysis.png']); %save figure

% set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
%     'PaperPosition', [0 0 N_column*13 N_row*12 ], 'PaperSize', [N_column*13 N_row*12 ] );
print(cur_fig, '-dpdf', [data.pathname data.filename(1:end-4) '_analysis.pdf']); %save figure

print(cur_fig, '-dpng', [data.pathname data.filename(1:end-4) '_analysis.png']); %save figure
disp('finished')

