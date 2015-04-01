%% startup
close all, clear all, clc

%% load gel data
gelData_raw = load_gel_image('data_dir', data_directory);

%%
[gelData_raw, cf] = check_gel_saturation(gelData_raw);

%% background correct data
gelData = background_correct_gel_image(gelData_raw, 'numberOfAreas', 4);

%% rotate image
gelData = rotate_gel_image(gelData);


%% overlay images 
[ch1_shift, ch1_dx, ch1_dy] = overlay_image(gelData.images{2}, gelData.images{1}, 'display', 'off');

gelData.images_raw = gelData.images;
gelData.images = {ch1_shift, gelData.images{2} };


%% integrate bands
bands = get_band_intensities(gelData);

%% create output dir
prefix_out = [gelData.filenames{1}(1:end-4) '_bands'];
prefix_out = [gelData.filenames{1}(1:end-4) '_bands-analysis_' datestr(now, 'yyyy-mm-dd_HH-MM')];
tmp = inputdlg({'Name of analysis (prefix):'}, 'Name of analysis (prefix):' , 1, {prefix_out} );
prefix_out = tmp{1};
path_out = [gelData.pathnames{1} prefix_out filesep];
mkdir(path_out);

%% save data
save([path_out prefix_out '_data.mat'])

%%
GFP_to_FS = zeros(size(bands.intensities,1), 2);
cd('/Users/jonasfunke/Documents/MATLAB/MATLAB_TOOLBOX/TYPHOON/private')
for i=1:size(bands.intensities,1)
    pos = bands.positions(i,:);
    subGFP = gelData.images{1}( pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3) );
    subFS = gelData.images{2}( pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3) );
    GFP_to_FS(i,:) = calculate_ration_of_areas(subGFP, subFS, 'display', 'off');
    
end
%%
close all
subplot(2,1,1)
plot(bands.intensities)
legend({'GFP', 'FS'})
subplot(2,1,2)
plot(1:size(bands.intensities,1), bands.intensities(:,1)./bands.intensities(:,2), '.-', ...
    1:size(bands.intensities,1), GFP_to_FS(:,1), '.-')
legend({'GFP/FS'})

%%
close all
j = [1:6; 7:12; 13:18; 19:24; 25:30];
hold all
for i=1:5
    %plot(1:6, bands.intensities(j(i,:),1)./bands.intensities(j(i,:),2), '.-', ...
    %1:6, GFP_to_FS(j(i,:),1), '.-')
    y = GFP_to_FS(j(i,:),1);
    
    %y = (y-y(1))./(y(5)-y(1));
    plot(1:6,y , '.-')
end
legend({'3', '3-26, 39 A', '3-198, 23 A', '3-212, 39 A', '132-212, 31 A'})

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

