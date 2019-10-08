%%
close all, clear all, clc

[filenames, pathname]=uigetfile('*.fcs','Select the fcs files','MultiSelect','on');

%%
for i=1:length(filenames)
    data(i) = load_fcs_data(pathname, filenames{i});

    
end

%%
fname = '2019-10-07-LGv7-2h-11.fcs';
pname = '/Users/jonasfunke/Dropbox (Personal)/Plectonic_Experiments/data_FACS/2019-10-07-LGv5 and LGv7/';
path_out = '/Users/jonasfunke/Dropbox (Personal)/Plectonic_Experiments/data_FACS/2019-10-07-LGv5 and LGv7';
prefix_out = '2019-10-07-LGv7-2h-11';

data = load_fcs_data(pname, fname);




%%



%% OLD STUFF
[fcsdat, fcshdr, fcsdatscaled, fcsdat_comp]  = fca_readfcs(fname);

%%
N_channel = fcshdr.NumOfPar;

i_plot = 2:2:16;

cur_fig = figure(1); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 20 10*length(i_plot) ], 'PaperSize', [20 10*length(i_plot)] );

for i=i_plot
    disp( [ num2str(i) ' - ' fcshdr.par(i).name])
    subplot(length(i_plot), 1, (i-2)/2+1)
    histogram(fcsdat(:,i))
    legend(fcshdr.par(i).name)
    title(['Median = ' num2str(median(fcsdat(:,i)))])
end

print(cur_fig, '-dpdf', [path_out filesep prefix_out '_raw.pdf']); %save figure


%% Select region of interest (ROI)
cur_fig = figure(2); clf

i=2; j=4; %FSC-A vs SSC-A
scatter(fcsdat(:,i), fcsdat(:,j), '.'), hold on
xlabel(fcshdr.par(i).name), ylabel(fcshdr.par(j).name)
set(gca,'xscale','log','yscale','log')
grid on
title('Select Region of Interest')

roi = drawpolygon;
i_gated = inROI(roi,fcsdat(:,i), fcsdat(:,j)) ;
scatter(fcsdat(i_gated,i), fcsdat(i_gated,j), 'r.') 

pause
close(2)

%% plot gated histograms
xlim = [0 1e4];
cur_fig = figure(1); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 20 10*length(i_plot) ], 'PaperSize', [20 10*length(i_plot)] );

for i=i_plot
    disp( [ num2str(i) ' - ' fcshdr.par(i).name])
    subplot(N_channel-1, 1, i)
    histogram(fcsdat(i_gated,i))
    legend(fcshdr.par(i).name)
    title(['Median = ' num2str(median(fcsdat(i_gated,i)))])
    %set(gca, 'XLim', xlim)
end

print(cur_fig, '-dpdf', [path_out filesep prefix_out '_gated.pdf']); %save figure


%% PLot scatterplot

cur_fig = figure(2); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 30 30 ], 'PaperSize', [30 30] );

i = 2;
for j=4:2:16
    disp( [ num2str(j) ' - ' fcshdr.par(j).name])
    subplot(3, 3, (j-4)/2+1)
    scatter(fcsdat(:,i), fcsdat(:,j),1, '.'), hold on
    scatter(fcsdat(i_gated,i), fcsdat(i_gated,j), 1,'.') 
    xlabel(fcshdr.par(i).name), ylabel(fcshdr.par(j).name)
    set(gca,'xscale','log','yscale','log')
    grid on
    legend({'Raw', 'Gated'}, 'location', 'NorthWest')
end
print(cur_fig, '-dpdf', [path_out filesep prefix_out '_scatter.pdf']); %save figure


%% PLot FL5-A channel histogram
i = 14; % FL5-A
cur_fig = figure(3); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 20 10 ], 'PaperSize', [20 10] );

histogram(fcsdat(i_gated,i))
title(['Median = ' num2str(median(fcsdat(i_gated,i)))])
xlabel(fcshdr.par(i).name), ylabel('Counts')
%set(gca, 'XLim', xlim)

print(cur_fig, '-dpdf', [path_out filesep prefix_out '_' fcshdr.par(i).name '_gated.pdf']); %save figure


