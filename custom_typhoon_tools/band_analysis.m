%% startup
close all, clear all, clc

%% load gel data
gelData_raw = load_gel_image('data_dir', data_directory);

%% background correct data
gelData = background_correct_gel_image(gelData_raw, 'numberOfAreas', 4);

%% rotate image
gelData = rotate_gel_image(gelData);

%% integrate bands
bands = get_band_intensities(gelData);

%% create output dir
prefix_out = [gelData.filenames{1}(1:end-4) '_bands'];
tmp = inputdlg({'Name of analysis (prefix):'}, 'Name of analysis (prefix):' , 1, {prefix_out} );
prefix_out = tmp{1};
path_out = [gelData.pathnames{1} prefix_out filesep];
mkdir(path_out);

%% save data
save([path_out prefix_out '_data.mat'])

%% compute yield
mono1 = bands.intensities(2:3:end,1);
dimer1 = bands.intensities(1:3:end,1);

mono2 = bands.intensities(3:3:end,2);
dimer2 = bands.intensities(1:3:end,2);
yield = [ [dimer1 ./ (dimer1 + mono1)],  [dimer2 ./ (dimer2 + mono2)]];

close all
plot(yield)

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

