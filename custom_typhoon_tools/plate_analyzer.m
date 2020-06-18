%% startup
close all, clear all, clc

%% load gel data
gelData_raw = load_gel_image('data_dir', data_directory, 'n_images', 3);

%% check for saturation
gelData_raw = check_gel_saturation(gelData_raw);

%%
gelData = gelData_raw;
%% overlay images 
[dd_shift, dd_dx, dd_dy] = overlay_image(gelData.images{2}, gelData.images{1}, 'display', 'off');
[da_shift, da_dx, da_dy] = overlay_image(gelData.images{2}, gelData.images{3}, 'display', 'off');

gelData.images_raw = gelData.images;
gelData.images = {dd_shift, gelData.images{2}, da_shift };

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
%%

dd = bandData.intensities(:,1);
aa = bandData.intensities(:,2);
da = bandData.intensities(:,3);

i_bg = [1,2];
dd = dd-mean(dd(i_bg));
aa = aa-mean(aa(i_bg));
da = da-mean(da(i_bg));

i_donly = 3;
i_aonly = 4;
leak = da(i_donly)/dd(i_donly);
dir = da(i_aonly)/aa(i_aonly);
%da = da - leak*dd - dir*aa;

close all
subplot(2,1,1)
plot(1:24, dd, 'g.-', ...
    1:24, aa, 'r.-', ...
    1:24, da, 'b.-')
grid on

subplot(2,1,2)
plot(1:24, da./(da+dd), 'k.-')
grid on



%% save data
close all
disp('Saving data...')
save([path_out prefix_out '_data.mat'])
disp('data saved...')


%% plot areas 
n_bands = size(bandData.positions,1);
areas =  bandData.positions;
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','points','PaperPosition', [0 0 1000 500], 'Position', [0 1000 1000 500]);
imagesc(gelData.images{2}, [0 3.*std(gelData.images{1}(:))]), axis image, colormap gray, hold on

for i=1:n_bands
    rectangle('Position', areas(i,:), 'EdgeColor', 'r', 'Linewidth', 1);
    text(areas(i,1)+areas(i,3)/2, areas(i,2) , num2str(i), 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 8)
end
set(gca, 'XTickLabel', [], 'YTickLabel', [])
print(cur_fig, '-dtiff', '-r 500' , [path_out filesep 'bands.tif']); %save figure

%% Plot 
n_bands = size(bandData.intensities,1);
cur_fig = figure;
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 20 25]);


subplot(5, 1, 1)
%plot(1:n_bands, bandData.intensities(:,2), 'r.-', 1:n_bands, bandData.intensities(:,1), 'g.-', 1:n_bands, bandData.intensities(:,4), 'b.-')
plot( 1:n_bands, gamma_calc.*bandData.intensities(:,1)./bandData.intensities(:,2), 'g.--', 1:n_bands, bandData.intensities(:,3)./bandData.intensities(:,2), 'b.--', ...
     1:n_bands, gamma_calc.*DD_div_AA(:,1), 'g.-', 1:n_bands, DA_div_AA(:,1), 'b.-')
xlabel('Lane'), ylabel('Relative band intensity')
legend({'int D->D/A->A', 'int D->A/A->A', 'scat D->D/A->A', 'scat D->A/A->A'}, 'location', 'best')
set(gca, 'XLim', [1 n_bands]);

subplot(5, 1, 2)
plot(1:n_bands, bandData.intensities(:,1)/max(bandData.intensities(:,1)), 'g.-', ...
    1:n_bands, bandData.intensities(:,2)/max(bandData.intensities(:,2)), 'r.-', ...+
    1:n_bands, bandData.intensities(:,3)/max(bandData.intensities(:,3)), 'b.-')
xlabel('Lane'), ylabel('Normalized band intensity')
legend({'D->D', 'A->A',  'D->A'}, 'location', 'west')
set(gca, 'XLim', [1 n_bands]);

subplot(5, 1, 3)
plot(1:n_bands, (bandData.intensities(:,1)./bandData.intensities(:,2))/max(bandData.intensities(:,1)./bandData.intensities(:,2)), 'g.-', ...
    1:n_bands, (bandData.intensities(:,2))/max(bandData.intensities(:,2)), 'r.-', ...+
    1:n_bands, (bandData.intensities(:,3)./bandData.intensities(:,2))/max(bandData.intensities(:,3)./bandData.intensities(:,2)), 'b.-')
xlabel('Lane'), ylabel('Normalized rel. band intensity')
legend({'D->D/A->A', 'A->A',  'D->A/A->A'}, 'location', 'west')
set(gca, 'XLim', [1 n_bands]);



subplot(5, 1, 4:5)
bar(1:n_bands, E), hold on
plot( 1:n_bands, E_integrate, 'k.--')
%plot( 1:n_bands, E_integrate, 'k.--', 1:n_bands, E, 'k.-')
xlabel('Lane'), ylabel('FRET efficiency')

title({['gamma=' num2str(gamma_calc) ]})
legend({  'FRET from scatterplot', 'FRET from intgration'})
set(gca, 'XLim', [1-0.5 n_bands+0.5], 'YLim', [0 1]);

print(cur_fig, '-dpdf', [path_out filesep 'FRET_normalized2.pdf']); %save figure

%% Plot 
close all
n_bands = size(bandData.intensities,1);
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 20 5], 'Position', [0 1000 2000 500]);

bar(  1:n_bands, E)
xlabel('Lane'), ylabel('FRET efficiency')

title({['gamma=' num2str(gamma_calc) ]})
set(gca, 'XLim', [0 n_bands+1], 'YLim', [0.0 1]);

print(cur_fig, '-dtiff', '-r 500' , [path_out filesep 'FRET_normalized_barplot.tif']); %save figure

print(cur_fig, '-depsc2' , [path_out filesep 'FRET_normalized_barplot.eps']); %save figure

