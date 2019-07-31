%% startup
close all, clear all, clc

%% load gel data
gelData_raw = load_gel_image('data_dir', data_directory, 'n_images', 2);

%%
[gelData_raw, cf] = check_gel_saturation(gelData_raw);

%% background correct data
gelData = background_correct_gel_image(gelData_raw, 'numberOfAreas', 4);

%% overlay images 
[ch1_shift, ch1_dx, ch1_dy] = overlay_image(gelData.images{2}, gelData.images{1}, 'display', 'off');

gelData.images_raw = gelData.images;
gelData.images = {ch1_shift, gelData.images{2} };


%% integrate bands
bands = get_band_intensities(gelData);

%% create output dir
prefix_out = [gelData.filenames{1}(1:end-4) '_bands-analysis_' datestr(now, 'yyyy-mm-dd_HH-MM')];
tmp = inputdlg({'Name of analysis (prefix):'}, 'Name of analysis (prefix):' , 1, {prefix_out} );
prefix_out = tmp{1};
path_out = [gelData.pathnames{1} prefix_out filesep];
mkdir(path_out);

%% save data
save([path_out prefix_out '_data.mat'])

%% plot areas 
areas =  bands.positions;
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeter','PaperPosition', [0 0 20 15 ], 'PaperSize', [20 15]);

subplot(2, 1, 1)
imagesc(gelData.images{1}, [0 3.*std(gelData.images{1}(:))]), axis image, colormap gray, hold on
for i=1:size(bands.intensities,1)
    rectangle('Position', areas(i,:), 'EdgeColor', 'r', 'Linewidth', 1);
    text(areas(i,1)+areas(i,3)/2, areas(i,2) , num2str(i), 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 8)
end
set(gca, 'XTickLabel', [], 'YTickLabel', [])

subplot(2, 1, 2)
imagesc(gelData.images{2}, [0 3.*std(gelData.images{2}(:))]), axis image, colormap gray, hold on
for i=1:size(bands.intensities,1)
    rectangle('Position', areas(i,:), 'EdgeColor', 'r', 'Linewidth', 1);
    text(areas(i,1)+areas(i,3)/2, areas(i,2) , num2str(i), 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 8)
end
set(gca, 'XTickLabel', [], 'YTickLabel', [])

print(cur_fig, '-dpdf' , [path_out filesep 'bands.pdf']); %save figure

%% Plot
n_bands = size(bands.intensities,1);




cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeter','PaperPosition', [0 0 15 30 ], 'PaperSize', [15 30]);

subplot(5, 1, 1)
imagesc(gelData.images{1}, [0 3.*std(gelData.images{1}(:))]), axis image, colormap gray, hold on
for i=1:size(bands.intensities,1)
    rectangle('Position', areas(i,:), 'EdgeColor', 'r', 'Linewidth', 1);
    text(areas(i,1)+areas(i,3)/2, areas(i,2) , num2str(i), 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 8)
end
set(gca, 'XTickLabel', [], 'YTickLabel', [])
title('Channel 1')

subplot(5, 1, 2)
imagesc(gelData.images{2}, [0 3.*std(gelData.images{2}(:))]), axis image, colormap gray, hold on
for i=1:size(bands.intensities,1)
    rectangle('Position', areas(i,:), 'EdgeColor', 'r', 'Linewidth', 1);
    text(areas(i,1)+areas(i,3)/2, areas(i,2) , num2str(i), 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 8)
end
set(gca, 'XTickLabel', [], 'YTickLabel', [])
title('Channel 2')

subplot(5, 1, 3)
plot(1:n_bands, bands.intensities(:,1), '.-', ...
    1:n_bands, bands.intensities(:,2), '.-' )
legend({'Channel 1', 'Channel 2'})
xlabel('Band')
ylabel('Mean intensity')
grid on

% calculate ratios
f1 = bands.intensities(:,2)./(bands.intensities(:,1)+bands.intensities(:,2));
r2 = zeros(size(f1));
for i=1:n_bands
    pos = bands.positions(i,:);
    band_ch1 = gelData.images{1}( pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3) );
    band_ch2 = gelData.images{2}( pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3) );
    [p_fit, ~, ~] = calculate_ration_of_areas(band_ch2, band_ch1, 'display', 'off');
    r2(i) = p_fit(1);
end
f2 = r2./(1+r2);


subplot(5, 1, 4)
plot(1:n_bands, f1, 'k.-', 1:n_bands, f2, 'k.--' )
legend({'integrated', 'scatter plot'})
ylabel('Fraction ch2/(ch1+ch2)')
xlabel('Band')
grid on

subplot(5, 1, 5)
plot(1:n_bands,f1/max(f1), 'k.-', 1:n_bands,f2/max(f2), 'k.--' )
legend({'integrated', 'scatter plot'})
ylabel('Normalized fraction ch2/(ch1+ch2) to max and min')
xlabel('Band')
grid on

print(cur_fig, '-dpdf' , [path_out filesep 'band_ratios.pdf']); %save figure


%%


%%
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeter','PaperPosition', [0 0 15 30 ], 'PaperSize', [15 30]);

subplot(3, 1, 1)
plot(1:n_bands, bands.intensities(:,1), '.-', ...
    1:n_bands, bands.intensities(:,2), '.-' )
legend({'Channel 1', 'Channel 2'})
xlabel('Band')
ylabel('Mean intensity')
grid on

% calculate ratios
f1 = bands.intensities(:,2)./(bands.intensities(:,1));
r2 = zeros(size(f1));
for i=1:n_bands
    pos = bands.positions(i,:);
    band_ch1 = gelData.images{1}( pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3) );
    band_ch2 = gelData.images{2}( pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3) );
    [p_fit, ~, ~] = calculate_ration_of_areas(band_ch2, band_ch1, 'display', 'off');
    r2(i) = p_fit(1);
end
f2 = r2;


subplot(3, 1, 2)
plot(1:n_bands, f1, 'k.-', 1:n_bands, f2, 'k.--' )
legend({'integrated', 'scatter plot'})
ylabel('Fraction ch2/ch1')
xlabel('Band')
grid on

subplot(3, 1, 3)
plot(1:n_bands,f1/max(f1), 'k.-', 1:n_bands,f2/max(f2), 'k.--' )
legend({'integrated', 'scatter plot'})
ylabel('Normalized fraction ch2/ch1 to max and min')
xlabel('Band')
grid on

print(cur_fig, '-dpdf' , [path_out filesep 'band_ratios_2.pdf']); %save figure
