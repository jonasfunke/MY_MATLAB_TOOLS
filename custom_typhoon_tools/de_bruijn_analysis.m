%% startup
close all, clear all, clc

%% load gel data
gelData_raw = load_gel_image('data_dir', data_directory, 'n_images', 2);

%% check for saturation
gelData_raw = check_gel_saturation(gelData_raw);

%% background correct data
gelData = background_correct_gel_image(gelData_raw, 'histogram_background', 'on');
gelData.background
gelData = background_correct_gel_image(gelData_raw, 'numberOfAreas', 4);
gelData.images_raw = gelData.images;


%% create output dir
prefix_out = [gelData.filenames{1}(1:end-4) '_analysis_' datestr(now, 'yyyy-mm-dd_HH-MM')];
tmp = inputdlg({'Name of analysis (prefix):'}, 'Name of analysis (prefix):' , 1, {prefix_out} );
prefix_out = tmp{1};
path_out = [gelData.pathnames{1} prefix_out filesep];
mkdir(path_out);

%% deterime profiles
profileData = get_gel_lanes(gelData, 'display', 'on', 'cutoff', 0.05, 'selection_type', 'automatic');
%profileData = get_gel_lanes(gelData, 'display', 'on', 'cutoff', 0.2, 'selection_type', 'manual');


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

myleg = {};
for i=1:length(profileData.profiles)
    plot(profileData.lanePositions(i,3):profileData.lanePositions(i,4), profileData.profiles{i}), hold on
    myleg = [myleg, {['Lane ' num2str(i)]}];
end
legend(myleg)
ylabel('Raw Intensity')
print(cur_fig, '-dtiff', '-r 300' , [path_out filesep 'Intensity.tif']); %save figure

%% Plot 
cur_fig = figure;

myleg = {};
for i=1:length(profileData.profiles)
    plot(profileData.lanePositions(i,3):profileData.lanePositions(i,4), profileData.profiles{i}/sum(profileData.profiles{i})), hold on
    myleg = [myleg, {['Lane ' num2str(i)]}];
end
legend(myleg)
ylabel('Normalized Intensity')

print(cur_fig, '-dtiff', '-r 300' , [path_out filesep 'Normalized_intensity.tif']); %save figure

%%

close all
max_values = zeros(length(profileData.profiles),2);


cur_fig = figure;
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 20 20]);
myleg = {};

for i=1:length(profileData.profiles)
    subplot(4,4,i)
    plot(profileData.lanePositions(i,3):profileData.lanePositions(i,4), profileData.profiles{1,i}, 'r', ...
        profileData.lanePositions(i,3):profileData.lanePositions(i,4), profileData.profiles{2,i}, 'g'), hold on    
    max_values(i,1) = max(profileData.profiles{1,i});
    max_values(i,2) = max(profileData.profiles{2,i});
    title(['Lane ' num2str(i)])
    ylabel('Raw Intensity')
xlabel('Migration distance [px]')
end

print(cur_fig, '-dpdf', [path_out filesep 'De-bruijn_analysis_profiles.pdf']); %save figure



%%
close all

cur_fig = figure;
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 20 5]);

plot(1:size(max_values, 1), max_values(:,1), '.-r', 1:size(max_values, 1), max_values(:,2), 'g.-')
set(gca, 'Xlim',[0.5 11.5], 'Xtick', [1:11], 'XtickLabel', myleg)
ylabel('Maximum intensity')
%xlabel('Lane')
legend({'de-bruijn', 'object'})



print(cur_fig, '-dpdf', [path_out filesep 'De-bruijn_analysis_intensities.pdf']); %save figure



%%

defect_rate = zeros(size(max_values, 1),2);
w_band = 15;
for i=1:length(profileData.profiles)
    [I_max, i_max] = max(profileData.profiles{2, i}); % Find maximum of lane based on 2nd channel (object)

    pos = profileData.lanePositions(i,:); 
    areas(i,:) = [pos(1), pos(3)+i_max-w_band, pos(2)-pos(1), 2*w_band];

    sub_debruijn = gelData.images{1}(pos(3)+i_max-w_band:pos(3)+i_max+w_band , pos(1):pos(2));
    sub_object = gelData.images{2}(pos(3)+i_max-w_band:pos(3)+i_max+w_band , pos(1):pos(2));
    
    defect_rate(i,:) = calculate_ration_of_areas(sub_debruijn, sub_object, 'display', 'off');
   

end

%%
cur_fig = figure;
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 20 10]);
plot( 1:size(max_values, 1), max_values(:,1)./max_values(:,2), 'k.-', 1:size(max_values, 1), defect_rate(:,1), 'k.--')
N_lanes = size(max_values, 1);
ylabel('Defect rate (de-bruijin/object) [a.u.]')
xlabel('Lane')
set(gca, 'Xlim',[0.5 N_lanes+0.5], 'Xtick', [1:N_lanes])
legend({'From maximum intenstiy', 'From scatter plot'}, 'location', 'northwest')

print(cur_fig, '-dpdf', [path_out filesep 'De-bruijn_analysis.pdf']); %save figure



%%

disp('done')