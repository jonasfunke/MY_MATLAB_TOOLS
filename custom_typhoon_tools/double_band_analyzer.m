%% startup
close all, clear all, clc

%% load gel data
gelData_raw = load_gel_image('data_dir', data_directory, 'n_images', 1);

%% check for saturation
gelData_raw = check_gel_saturation(gelData_raw);

%% background correct data
gelData = background_correct_gel_image(gelData_raw, 'numberOfAreas', 4);

%% create output dir
prefix_out = [gelData.filenames{1}(1:end-4) '_analysis_' datestr(now, 'yyyy-mm-dd_HH-MM')];
tmp = inputdlg({'Name of analysis (prefix):'}, 'Name of analysis (prefix):' , 1, {prefix_out} );
prefix_out = tmp{1};
path_out = [gelData.pathnames{1} prefix_out filesep];
mkdir(path_out);

%% integrate bands
bandData = get_band_intensities(gelData);
pause(0.1)
close all

%% fit double bands
% compute band profiles
bandData.profiles = cell(1, size(bandData.positions,1));
bandData.migration_distance = cell(1, size(bandData.positions,1));
for i=1:size(bandData.positions,1)
    tmp = mean(gelData.images{1}( bandData.positions(i,2):bandData.positions(i,2)+bandData.positions(i,4), bandData.positions(i,1):bandData.positions(i,1)+bandData.positions(i,3) ), 2);
    bandData.profiles{i} = tmp;
    bandData.migration_distance{i} = double(bandData.positions(i,2):bandData.positions(i,2)+bandData.positions(i,4))';
end

%% fit band profile
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters',...
    'PaperPosition', [0 0 15 size(bandData.positions,1)*5 ], 'PaperSize', [15 size(bandData.positions,1)*5]);

bandData.fits = cell(1, size(bandData.positions,1));
bandData.fits_2 = cell(1, size(bandData.positions,1));
bandData.peaks = cell(2, size(bandData.positions,1));

ft_gauss1 = fittype('gauss1');
ft_gauss2_same_width = fittype('a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c1)^2)');
for i=1:size(bandData.positions,1)

    subplot(size(bandData.positions,1), 1, i), cla
    [~, index] = findpeaks(bandData.profiles{i}, 'SortStr','descend');
    if length(index)>1
        bandData.peaks{i} = sort(index(1:2));
    else
        bandData.peaks{i} = [index index];
    end
    
    
%     bandData.fits2{i} = fit(bandData.migration_distance{i}, bandData.profiles{i}, 'gauss2');
%     bandData.fits{i} = fit(bandData.migration_distance{i}, bandData.profiles{i}, 'gauss2', ...
%         'StartPoint', [bandData.profiles{i}(bandData.peaks{i}(1)), bandData.migration_distance{i}(bandData.peaks{i}(1)) , 10, ...
%         bandData.profiles{i}(bandData.peaks{i}(2)), bandData.migration_distance{i}(bandData.peaks{i}(2)) , 10 ]);
    bandData.fits3{i} = fit(bandData.migration_distance{i}, bandData.profiles{i}, ft_gauss2_same_width, ...
        'StartPoint', [bandData.profiles{i}(bandData.peaks{i}(1)), ...
        bandData.profiles{i}(bandData.peaks{i}(2)), bandData.migration_distance{i}(bandData.peaks{i}(1)) , bandData.migration_distance{i}(bandData.peaks{i}(2)), 10   ]);

    %bandData.fits{i}
    plot(bandData.migration_distance{i}, bandData.profiles{i}), hold on

%     plot(bandData.migration_distance{i}, bandData.profiles{i}, ...
%         bandData.migration_distance{i},  bandData.fits{i}(bandData.migration_distance{i}), '--', ...
%         bandData.migration_distance{i},  bandData.fits2{i}(bandData.migration_distance{i}), '--')
    plot(bandData.migration_distance{i}, bandData.fits3{i}(bandData.migration_distance{i}), ':', ...
        bandData.migration_distance{i},  ft_gauss1(bandData.fits3{i}.a1, bandData.fits3{i}.b1, bandData.fits3{i}.c1,bandData.migration_distance{i}) , '-', ...
        bandData.migration_distance{i},  ft_gauss1(bandData.fits3{i}.a2, bandData.fits3{i}.b2, bandData.fits3{i}.c1,bandData.migration_distance{i}) , '-')
    plot(bandData.migration_distance{i}(bandData.peaks{i}) , bandData.profiles{i}(bandData.peaks{i}), 'k.'), hold on

  %  legend({'profile', ...
  %      ['fit w1=' num2str(round(bandData.fits{i}.c1)) ' w2=' num2str(round(bandData.fits{i}.c2)) ],...
  %      'peaks' })
  
  
    
end
print(cur_fig, '-dpdf' , [path_out filesep 'Profiles.pdf']); %save figure


%% save data
close all
disp('Saving data...')
save([path_out prefix_out '_data.mat'])
disp('data saved...')


%% plot areas 
n_bands = size(bandData.positions,1);
areas =  bandData.positions;
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','points','PaperPosition', [0 0 1000 500], 'Position', [0 1000 1000 500]);
imagesc(gelData.images{1}, [0 3.*std(gelData.images{1}(:))]), axis image, colormap gray, hold on

for i=1:n_bands
    rectangle('Position', areas(i,:), 'EdgeColor', 'r', 'Linewidth', 1);
    text(areas(i,1)+areas(i,3)/2, areas(i,2) , num2str(i), 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 8)
    
    %plot((areas(i,1)+areas(i,3)/2)*[1;1], bandData.migration_distance{i}(bandData.peaks{i}), 'g.')
    
    % plot alos positions from fit3
    plot((areas(i,1)+areas(i,3)/2)*[1;1], bandData.fits3{i}.b1, 'r.')
    plot((areas(i,1)+areas(i,3)/2)*[1;1], bandData.fits3{i}.b2, 'r.')

end
set(gca, 'XTickLabel', [], 'YTickLabel', [])
print(cur_fig, '-dtiff', '-r 500' , [path_out filesep 'bands.tif']); %save figure

%% Plot 
n_bands = size(bandData.intensities,1);



cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 20 15], 'PaperSize', [20 15]);

heights = zeros(n_bands, 2);
for i=1:n_bands
    h1 = bandData.profiles{i}(bandData.peaks{i}(1));
    h1 = bandData.fits3{i}.a1;
    h2 = bandData.profiles{i}(bandData.peaks{i}(2));
    h2 = bandData.fits3{i}.a2;
    if bandData.fits3{i}.b2 < bandData.fits3{i}.b1 % switch
        heights(i,:) = [h2 h1];
    else
        heights(i,:) = [h1 h2];
    end
end
subplot(2, 1, 1)
plot(1:n_bands, heights, '.-')
set(gca, 'XTick', 1:n_bands, 'Xlim', [0.5 n_bands+0.5])
grid on
legend({'slow', 'fast'}, 'location', 'best')
ylabel('Peak height')
xlabel('Double-band')

subplot(2, 1, 2)
%plot(1:n_bands, heights(:,1)./sum(heights,2), '.-', 1:n_bands, heights(:,2)./sum(heights,2), '.-')
bar(1:n_bands, heights(:,1)./sum(heights,2))
set(gca, 'XTick', 1:n_bands, 'Xlim', [0.5 n_bands+0.5], 'YLim', [-0.1 1.1])
%legend({'slow', 'fast'}, 'location', 'best')

grid on
ylabel('Fraction')
xlabel('Double-band')

set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 30 20 ], 'PaperSize', [30 20 ] );
print(cur_fig, '-dpdf', [path_out filesep 'Double-band_ratios.pdf']); %save figure

%%


subplot(2, 1, 2)
%plot(1:n_bands, heights(:,1)./sum(heights,2), '.-', 1:n_bands, heights(:,2)./sum(heights,2), '.-')
bar(1:n_bands, heights(:,1)./sum(heights,2))
set(gca, 'XTick', 1:n_bands, 'Xlim', [0.5 n_bands+0.5], 'YLim', [-0.1 1.1])
%legend({'slow', 'fast'}, 'location', 'best')

grid on
ylabel('Fraction')
xlabel('Double-band')

set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 30 20 ], 'PaperSize', [30 20 ] );
print(cur_fig, '-dpdf', [path_out filesep 'Double-band_ratios.pdf']); %save figure


%%

cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 20 7], 'PaperSize', [20 7]);



plot(1:n_bands, sum(heights,2)/mean(sum(heights,2)), '.-')
set(gca, 'XTick', 1:n_bands, 'Xlim', [0.5 n_bands+0.5], 'YLim', [0 1.5])

grid on
ylabel('Sum oh heights')
xlabel('Double-band')

print(cur_fig, '-dpdf', [path_out filesep 'Double-band_total.pdf']); %save figure
