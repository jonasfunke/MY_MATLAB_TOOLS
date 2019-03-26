%% startup
close all, clear all, clc

%% load gel data
gelData_raw = load_gel_image('data_dir', data_directory, 'n_images', 1);

%% check for saturation
gelData_raw = check_gel_saturation(gelData_raw);

%% background correct data
%gelData = background_correct_gel_image(gelData_raw, 'histogram_background', 'on');
%gelData.background
gelData = background_correct_gel_image(gelData_raw, 'numberOfAreas', 4);
gelData.images_raw = gelData.images;


%% create output dir
prefix_out = [gelData.filenames{1}(1:end-4) '_analysis_' datestr(now, 'yyyy-mm-dd_HH-MM')];
tmp = inputdlg({'Name of analysis (prefix):'}, 'Name of analysis (prefix):' , 1, {prefix_out} );
prefix_out = tmp{1};
path_out = [gelData.pathnames{1} prefix_out filesep];
mkdir(path_out);

%% deterime profiles
%profileData = get_gel_lanes(gelData, 'display', 'on', 'cutoff', 0.05, 'selection_type', 'automatic');
profileData = get_gel_lanes(gelData, 'display', 'on', 'cutoff', 0.1, 'selection_type', 'manual');


%% save data
close all
disp('Saving data...')
save([path_out prefix_out '_data.mat'])
disp('data saved...')

%% plot areas 
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','points','PaperPosition', [0 0 1000 500], 'Position', [0 1000 1000 500]);imagesc(gelData.images{1}, [0 3.*std(gelData.images{1}(:))]), axis image, colormap gray, hold on

areas = [profileData.lanePositions(:,1)   profileData.lanePositions(:,3) profileData.lanePositions(:,2)-profileData.lanePositions(:,1)  profileData.lanePositions(:,4)-profileData.lanePositions(:,3)];
for i=1:length(profileData.profiles)
    rectangle('Position', areas(i,:), 'EdgeColor', 'r', 'Linewidth', 1);
    text(areas(i,1)+areas(i,3)/2, areas(i,2) , num2str(i), 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 8)
end
set(gca, 'XTickLabel', [], 'YTickLabel', [])
print(cur_fig, '-dtiff', '-r 500' , [path_out filesep 'Lanes.tif']); %save figure


%% Plot 
cur_fig = figure;
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 20 10]);

myleg = {};
for i=1:length(profileData.profiles)
    plot(profileData.lanePositions(i,3):profileData.lanePositions(i,4), profileData.profiles{i}), hold on
    myleg = [myleg, {['Lane ' num2str(i)]}];
end
legend(myleg)
ylabel('Raw Intensity')
xlabel('Migration distance [px]')
print(cur_fig, '-dpdf', [path_out filesep 'Intensity.pdf']); %save figure

%% Plot 
cur_fig = figure;

myleg = {};
for i=1:length(profileData.profiles)
    plot(profileData.lanePositions(i,3):profileData.lanePositions(i,4), profileData.profiles{i}/sum(profileData.profiles{i})), hold on
    myleg = [myleg, {['Lane ' num2str(i)]}];
end
legend(myleg)
ylabel('Normalized Intensity')
xlabel('Migration distance [px]')
print(cur_fig, '-dtiff', '-r 300' , [path_out filesep 'Normalized_intensity.tif']); %save figure

%% Plot 
cur_fig = figure;


myleg = {};
for i=1:length(profileData.profiles)
    myleg = [myleg, {['Lane ' num2str(i)]}];
end

for i=1:length(profileData.profiles)
    [max_intensity, max_index] = max(profileData.profiles{i}/sum(profileData.profiles{i}));
    plot([profileData.lanePositions(i,3):profileData.lanePositions(i,4)]-max_index, profileData.profiles{i}/sum(profileData.profiles{i})/max_intensity), hold on
    %pause
end
legend(myleg)
    

legend(myleg(1:length(profileData.profiles)))
ylabel('Normalized Intensity')

print(cur_fig, '-dtiff', '-r 300' , [path_out filesep 'Normalized_intensity_to_max_shifted.tif']); %save figure



%%

cur_fig = figure;

for i=1:length(profileData.profiles)
    plot(profileData.lanePositions(i,3):profileData.lanePositions(i,4), profileData.profiles{i}/sum(profileData.profiles{i})), hold on
end
legend(myleg([ 3, 8:11]))
ylabel('Normalized Intensity')


disp('done')