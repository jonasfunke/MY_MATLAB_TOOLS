%% startup
close all, clear all, clc

bands = get_band_intensities('data_dir', data_directory, 'areas_bg', 4);

%%
save([bands.pathnames{1} 'band_data.mat'])
path_out = bands.pathnames{1} ;

%% plot
close all
cur_fig = figure();

subplot(2, 1, 1)
plot([10:51],bands.band_intensities(:,1), 'g', [10:51],bands.band_intensities(:,2), 'r'  )
legend({'GFP-channel', 'FS-channel'})
xlabel('Length of spacer [bp]')
ylabel('Intensity')

subplot(2, 1, 2)
plot([10:51], bands.band_intensities(:,1) ./ bands.band_intensities(:,2), 'k'  )
xlabel('Length of spacer [bp]')
ylabel('relative intensity(GFP/FS)')

print(cur_fig, '-dtiff','-r500' , [path_out filesep 'Raw_intensities.tif']); %save figure


%% plot normalized
close all
cur_fig = figure();

y1 = bands.band_intensities(:,1) ./ mean(bands.band_intensities(:,1));
y2 = bands.band_intensities(:,2) ./ mean(bands.band_intensities(:,2));

subplot(2, 1, 1)
plot([10:51],y1, 'g', [10:51],y2, 'r'  )
legend({'GFP-channel', 'FS-channel'})
xlabel('Length of spacer [bp]')
ylabel('Normalized intensity to mean')

subplot(2, 1, 2)
plot([10:51], y1 ./ y2, 'k'  )
xlabel('Length of spacer [bp]')
ylabel('relative intensity(GFP/FS)')

print(cur_fig, '-dtiff','-r500' , [path_out filesep 'normalized_intensities.tif']); %save figure

