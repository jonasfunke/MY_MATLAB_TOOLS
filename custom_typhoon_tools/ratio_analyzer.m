%% startup
close all, clear all, clc

%% load gel data
gelData_raw = load_gel_image('data_dir', data_directory, 'n_images', 2);

%% check for saturation
gelData_raw = check_gel_saturation(gelData_raw);

%% background correct data
gelData = background_correct_gel_image(gelData_raw, 'numberOfAreas', 4);

%% overlay images 
[ch2_shift, ch2_dx, ch2_dy] = overlay_image(gelData.images{1}, gelData.images{2}, 'display', 'off');
gelData.images_raw = gelData.images;
gelData.images = {gelData.images{1}, ch2_shift};

%% create output dir
prefix_out = [gelData.filenames{1}(1:end-4) '_analysis_' datestr(now, 'yyyy-mm-dd_HH-MM')];
tmp = inputdlg({'Name of analysis (prefix):'}, 'Name of analysis (prefix):' , 1, {prefix_out} );
prefix_out = tmp{1};
path_out = [gelData.pathnames{1} prefix_out filesep];
mkdir(path_out);

%% deterime profiles
profileData = get_gel_lanes(gelData, 'display', 'off', 'cutoff', 0.01);

%% Calculate ratio of bands
w_band = 10;
n_bands = size(profileData.lanePositions, 1);
I_mean = zeros(n_bands, gelData.nrImages);
ratio = zeros(n_bands, 1);

path0 = cd;
cd('/Users/jonasfunke/Documents/MATLAB/MATLAB_TOOLBOX/TYPHOON/private')
for i=1:n_bands
    [I_max, i_max] = max(profileData.profiles{1, i}); % Find maximum of lane based on 1st channel
    pos = profileData.lanePositions(i,:);
    
    band_ch1 = gelData.images{1}(pos(3)+i_max-w_band:pos(3)+i_max+w_band , pos(1):pos(2));
    band_ch2 = gelData.images{2}(pos(3)+i_max-w_band:pos(3)+i_max+w_band , pos(1):pos(2));
    p = calculate_ration_of_areas(band_ch2, band_ch1, 'display', 'off');
    ratio(i) = p(1);
    
    for channel = 1:gelData.nrImages
        cur_band =  gelData.images{channel}(pos(3)+i_max-w_band:pos(3)+i_max+w_band , pos(1):pos(2));
        I_mean(i,channel) = mean(cur_band(:));
    end
end
cd(path0)


%% plot areas 
%{
close all
imagesc(gelData.images{1}), axis image, colormap gray, hold on

for i=1:n_bands
    [I_max, i_max] = max(profileData.profiles{ 1, i}); % based on 1st channel
    pos = profileData.lanePositions(i,:);
    rectangle('Position', [pos(1), pos(3)+i_max-w, pos(2)-pos(1), 2*w], 'EdgeColor', 'r');
end
%}

%% save data
close all
disp('Saving data...')
save([path_out prefix_out '_data.mat'])
disp('data saved...')

%% Plot mean intensity of band
cur_fig = figure();
subplot(2,1,1)
plot(1:n_bands, I_mean(:,1), '.-', 1:n_bands,I_mean(:,2), '.-')
set(gca, 'XLim', [0, n_bands+1])
xlabel('Lane')
ylabel('Mean intensity [a.u.]')
legend({'channel 1', 'channel 2'}, 'location', 'best')


subplot(2,1,2)
plot(1:n_bands, I_mean(:,1)-mean(I_mean(:,1)), '.-', 1:n_bands,I_mean(:,2)-mean(I_mean(:,2)), '.-')
set(gca, 'XLim', [0, n_bands+1])
xlabel('Lane')
ylabel('Mean intensity - mean intensity of bands  [a.u.]')
legend({'channel 1', 'channel 2'}, 'location', 'best')

print(cur_fig, '-dtiff', '-r 500' , [path_out filesep 'Total_intensity.tif']); %save figure

%% plot Ratio
cur_fig = figure();
plot(1:n_bands, I_mean(:,2)./I_mean(:,1), '.-', 1:n_bands, ratio, '.-')
set(gca, 'XLim', [0, n_bands+1])
xlabel('Lane')
ylabel('Relative intensity [channel 2 / channel 1]')
legend({'From mean intensities', 'from scatter-plot'}, 'location', 'best')

print(cur_fig, '-dtiff', '-r 500' , [path_out filesep 'Relative_intensity.tif']); %save figure

%%
disp('done')