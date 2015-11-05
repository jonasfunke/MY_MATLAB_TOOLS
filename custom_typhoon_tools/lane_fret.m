%% startup
close all, clear all, clc

%% load gel data
gelData_raw = load_gel_image('data_dir', data_directory, 'n_images', 3);

%% check for saturation
gelData_raw = check_gel_saturation(gelData_raw);

%% background correct data
gelData = background_correct_gel_image(gelData_raw, 'histogram_background', 'on');
gelData.background
gelData = background_correct_gel_image(gelData_raw, 'numberOfAreas', 4);

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

%% leakage and direct-excitation correction factors
[leak_dir, leakdir_fig]  = calculate_leakage_directexcitation(gelData.images{1}, gelData.images{3}, gelData.images{2}, 'display', 'on');
print(leakdir_fig, '-dtiff', '-r 500' , [path_out filesep 'Leakage_DirEx_correction.tif']); %save figure

% write leakdir
save([path_out filesep 'Leakage_DirEx_data.txt'], 'leak_dir' ,'-ascii')

da_cor = gelData.images{3} - leak_dir(1,1).*gelData.images{1} - leak_dir(2,1).*gelData.images{2};%-leak_dir(1,2)-leak_dir(2,2); %correct image
gelData.images{3} = da_cor; % append to images

%% deterime profiles
profileData = get_gel_lanes(gelData, 'display', 'on', 'cutoff', 0.05, 'selection_type', 'automatic');

%% Calculate FRET efficiency for leading band
w_band = 10;
n_bands = size(profileData.lanePositions, 1);
I_mean = zeros(n_bands, gelData.nrImages);

areas = zeros(n_bands, 4);
DD_div_DA = zeros(n_bands, 2);
DD_div_AA = zeros(n_bands, 2);
DA_div_AA = zeros(n_bands, 2);

subplot(3,1,1)
imagesc(gelData.images{2}), axis image, colormap gray
hold on

path0 = cd;
for i=1:n_bands
    [I_max, i_max] = max(profileData.profiles{3, i}); % Find maximum of lane based on 2nd channel (D->A)
    pos = profileData.lanePositions(i,:);
    
    subplot(3,1,1)
    rectangle('Position', [pos(1), pos(3)+i_max-w_band, pos(2)-pos(1), 2*w_band], 'EdgeColor', 'r');
    
    areas(i,:) = [pos(1), pos(3)+i_max-w_band, pos(2)-pos(1), 2*w_band];
    subplot(3,1,2:3)
    hold off
    plot(1:length( profileData.profiles{1,i}), profileData.profiles{1,i}, 'g'), hold on
    plot(1:length( profileData.profiles{2,i}), profileData.profiles{2,i}, 'r'), hold on
    plot(1:length( profileData.profiles{3,i}), profileData.profiles{3,i}, 'b'), hold on
    vline(i_max, {'k'});
    vline(i_max-w_band, {'k--'});
    vline(i_max+w_band, {'k--'});
    
    %pause
    subDD = gelData.images{1}(pos(3)+i_max-w_band:pos(3)+i_max+w_band , pos(1):pos(2));
    subAA = gelData.images{2}(pos(3)+i_max-w_band:pos(3)+i_max+w_band , pos(1):pos(2));
    subDA = gelData.images{3}(pos(3)+i_max-w_band:pos(3)+i_max+w_band , pos(1):pos(2));
    DD_div_DA(i,:) = calculate_ration_of_areas(subDD, subDA, 'display', 'off');
   % pause
   % close all
    DD_div_AA(i,:) = calculate_ration_of_areas(subDD, subAA, 'display', 'off');
    DA_div_AA(i,:) = calculate_ration_of_areas(subDA, subAA, 'display', 'off');
    
    %pause
    %p = calculate_ration_of_areas(band_ch2, band_ch1, 'display', 'off');
    %ratio(i) = p(1);
    
    for channel = 1:gelData.nrImages
        cur_band =  gelData.images{channel}(pos(3)+i_max-w_band:pos(3)+i_max+w_band , pos(1):pos(2));
        I_mean(i,channel) = mean(cur_band(:));
    end
end
cd(path0)

close all

%%
close all
imagesc(gelData.images{2}), colormap gray, axis image

for i=1:n_bands
    rectangle('Position', areas(i,:), 'EdgeColor', 'r', 'Linewidth', 1);
    text(areas(i,1)+areas(i,3)/2, areas(i,2) , num2str(i), 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 8)
end
%%
x = inputdlg('Enter space-separated index for reference samples:',...
             'Sample', [1 50]);
i_gamma = str2num(x{:})
%%

% %i_gamma = [9];
% 
 E_soll = 0.5;
% %gamma_calc =  bandData.intensities(i_gamma,4).*(1./0.5 - 1) ./  bandData.intensities(i_gamma,1) 
%     
 E_raw = 1./(1+DD_div_DA(:,1));
 gamma_calc = (1-E_soll)./E_soll./DD_div_DA(i_gamma,1)
 gamma_calc = mean(gamma_calc)
%gamma_calc = 1;
E = 1./(1+gamma_calc.*DD_div_DA(:,1));

%%
%gamma_calc_integrate =  I_mean(i_gamma,3).*(1./E_soll - 1) ./  I_mean(i_gamma,1) 
E_integrate = I_mean(:,3) ./ (gamma_calc.*I_mean(:,1) + I_mean(:,3));



%% save data
close all
disp('Saving data...')
save([path_out prefix_out '_data.mat'])
disp('data saved...')


%% write corrected images
disp('Writing images')

t = Tiff([path_out filesep 'da_cor.tif'],'w');
t.setTag('Photometric',Tiff.Photometric.MinIsWhite);
t.setTag('BitsPerSample',16);
t.setTag('SampleFormat',Tiff.SampleFormat.UInt);
t.setTag('ImageLength',size(da_cor,1));
t.setTag('ImageWidth',size(da_cor,2));
t.setTag('SamplesPerPixel',1);
t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
t.write( uint16(da_cor-min(da_cor(:)))  );
t.close();

t = Tiff([path_out filesep 'da_cor+bg.tif'],'w');
t.setTag('Photometric',Tiff.Photometric.MinIsWhite);
t.setTag('BitsPerSample',16);
t.setTag('SampleFormat',Tiff.SampleFormat.UInt);
t.setTag('ImageLength',size(da_cor,1));
t.setTag('ImageWidth',size(da_cor,2));
t.setTag('SamplesPerPixel',1);
t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
t.write( uint16(da_cor+gelData.background{3}.p00)  );
t.close();
 
%% plot areas 

cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','points','PaperPosition', [0 0 1000 500], 'Position', [0 1000 1000 500]);imagesc(gelData.images{1}, [0 3.*std(gelData.images{1}(:))]), axis image, colormap gray, hold on

for i=1:n_bands
    rectangle('Position', areas(i,:), 'EdgeColor', 'r', 'Linewidth', 1);
    text(areas(i,1)+areas(i,3)/2, areas(i,2) , num2str(i), 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 8)
end
set(gca, 'XTickLabel', [], 'YTickLabel', [])
print(cur_fig, '-dtiff', '-r 500' , [path_out filesep 'bands.tif']); %save figure


%% Plot 
cur_fig = figure;
subplot(3, 1, 1)
%plot(1:n_bands, bandData.intensities(:,2), 'r.-', 1:n_bands, bandData.intensities(:,1), 'g.-', 1:n_bands, bandData.intensities(:,4), 'b.-')
plot( 1:n_bands, gamma_calc.*I_mean(:,1)./I_mean(:,2), 'g.--', 1:n_bands, I_mean(:,3)./I_mean(:,2), 'b.--', ...
     1:n_bands, gamma_calc.*DD_div_AA(:,1), 'g.-', 1:n_bands, DA_div_AA(:,1), 'b.-')
xlabel('Lane'), ylabel('Normalized bandintensity')
legend({'gamma * D->D / A->A', 'D->A / A->A'}, 'location', 'best')
set(gca, 'XLim', [1 n_bands],  'YLim', [0.1 0.4]);

subplot(3, 1, 2:3)
plot( 1:n_bands, E_integrate, 'k.--', 1:n_bands, E, 'k.-')
xlabel('Lane'), ylabel('FRET efficiency')

title({['gamma=' num2str(gamma_calc) ]})
legend({ 'FRET from intgration', 'FRET from scatterplot'})
set(gca, 'XLim', [1 n_bands], 'YLim', [0.3 0.6]);

print(cur_fig, '-dtiff', '-r 500' , [path_out filesep 'FRET_normalized.tif']); %save figure

%% Plot 
close all
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 20 5], 'Position', [0 1000 2000 500]);

bar(  1:n_bands, E)
xlabel('Lane'), ylabel('FRET efficiency')

title({['gamma=' num2str(gamma_calc) ]})
set(gca, 'XLim', [0 n_bands+1], 'YLim', [0.2 0.6]);

print(cur_fig, '-dtiff', '-r 500' , [path_out filesep 'FRET_normalized_barplot.tif']); %save figure

print(cur_fig, '-depsc2' , [path_out filesep 'FRET_normalized_barplot.eps']); %save figure

%%
disp('done')