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

save([path_out prefix_out '_data.mat'], 'yield' ,'-append')

%%
n_bands = size(bands.intensities,1);
n_samples = n_bands / 3;


%% plot total intensity
close all
cur_fig = figure()
plot(1:n_samples, mono1+dimer1, 'g.-', 1:n_samples, mono2+dimer2, 'r.-')
xlabel('Lane')
ylabel('Total Intensity')
legend({'cy3-channel', 'cy5-channel'}, 'Location', 'best')
set(gca, 'XLim', [1 n_samples])
print(cur_fig, '-dtiff', '-r 500' , [path_out filesep prefix_out '_total-intensity.tif']); %save figure

%% plot yield
close all
cur_fig = figure();
plot(1:n_samples, yield(:,1), 'g.-', 1:n_samples, yield(:,2), 'r.-')
set(gca, 'YLim', [0 0.5])
xlabel('Lane')
ylabel('Yield')
legend({'cy3-channel', 'cy5-channel'}, 'Location', 'best')
set(gca, 'XLim', [1 n_samples])

print(cur_fig, '-dtiff', '-r 500' , [path_out filesep prefix_out '_yield.tif']); %save figure


%%

subplot(1,2,1)
plot(mono1+dimer1, yield(:,1), 'g.')

subplot(1,2,2)
plot(mono2+dimer2, yield(:,2), 'r.')

%%
disp('Finished')
