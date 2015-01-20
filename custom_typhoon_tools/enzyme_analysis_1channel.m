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

%%
n_bands = size(bands.intensities,1);
n_samples = n_bands / 2;

%% compute yield
mono = bands.intensities(2:2:end,1);
dimer = bands.intensities(1:2:end,1);

yield =  [dimer ./ (dimer + mono)];

%% save data
save([path_out prefix_out '_data.mat'])



%% plot total intensity
close all
cur_fig = figure()
plot( 1:n_samples, mono+dimer, 'r.-')
xlabel('Lane')
ylabel('Total Intensity')
legend({'cy5-channel'}, 'Location', 'best')
set(gca, 'XLim', [1 n_samples], 'YLim', [0 1.2*max(mono+dimer)])
print(cur_fig, '-dtiff', '-r 500' , [path_out filesep prefix_out '_total-intensity.tif']); %save figure

%% plot yield
cur_fig = figure();
plot(1:n_samples, yield, 'r.-')
set(gca, 'YLim', [0 0.5])
xlabel('Lane')
ylabel('Yield')
legend({'cy3-channel', 'cy5-channel'}, 'Location', 'best')
set(gca, 'XLim', [1 n_samples])

print(cur_fig, '-dtiff', '-r 500' , [path_out filesep prefix_out '_yield.tif']); %save figure


%%
cur_fig = figure();

cur_fig = figure()
subplot(2,1,1)
plot( 1:n_samples, mono+dimer, 'r.-')
xlabel('Lane')
ylabel('Total Intensity')
legend({'cy5-channel'}, 'Location', 'best')
set(gca, 'XLim', [1 n_samples], 'YLim', [0 1.2*max(mono+dimer)])

subplot(2,1,2)
plot(1:n_samples, yield, 'r.-')
set(gca, 'YLim', [0 0.5])
xlabel('Lane')
ylabel('Yield')
legend({'cy3-channel', 'cy5-channel'}, 'Location', 'best')
set(gca, 'XLim', [1 n_samples])

%%
%close all
%plot(mono+dimer, yield, 'r.')
%R = corrcoef(mono+dimer, yield)
%%
disp('Finished')
