%% startup
close all, clear all, clc

%% load gel data
gelData_raw = load_gel_image('data_dir', data_directory, 'n_images', 3);

%% check for saturation
gelData_raw = check_gel_saturation(gelData_raw);

%% create output dir
prefix_out = [gelData_raw.filenames{1}(1:end-4) '_analysis_' datestr(now, 'yyyy-mm-dd_HH-MM')];
tmp = inputdlg({'Name of analysis (prefix):'}, 'Name of analysis (prefix):' , 1, {prefix_out} );
prefix_out = tmp{1};
path_out = [gelData_raw.pathnames{1} prefix_out filesep];
mkdir(path_out);
%%
N_sample = 36;

[tmp, bandPositions] = get_area_intensities(gelData_raw.images, N_sample, 'resizable', false); %cell of images, number of bands, 

%%
%close all
N_sample = size(tmp,3);
dd = zeros(N_sample,1);
aa = zeros(N_sample,1);
da = zeros(N_sample,1);
for i=1:size(tmp,3)
    tmp2 = tmp(:,:,i,1);
    dd(i) = median(tmp2(:));
    
    tmp3 = tmp(:,:,i,2);
    aa(i) = median(tmp3(:));
    
    tmp4 = tmp(:,:,i,3);
    da(i) = median(tmp4(:));
    
    
%     subplot(1, 3, 1), cla
%     hist(tmp2(:)), hold on
%     vline(mean(tmp2(:)), {'r--'});
%     vline(median(tmp2(:)), {'r'});
%     
%     subplot(1, 3, 2), cla
%     hist(tmp3(:)), hold on
%     vline(mean(tmp3(:)), {'r--'});
%     vline(median(tmp3(:)), {'r'});
%     
%     subplot(1, 3, 3), cla
%     hist(tmp4(:)), hold on
%     vline(mean(tmp4(:)), {'r--'});
%     vline(median(tmp4(:)), {'r'});
%     pause
end




%

% correct backgroud
i_bg = [1 2 9 10 17 18 19 20];
titles ={ 'FoB6','FoB6','D 1nM','A 1nM','D 5nM','A 5nM','D 10nM','A 10nM', ...
    'FoB6','FoB6','LF 1nM','HF 1nM','LF 5nM','HF 5nM','LF 10nM','HF 10nM', ...
    'PBS5', 'PBS5', 'PBS5', 'PBS5', 'C', 'C', 'C', 'C', ...
    'D 1nM+C','A 1nM+C','D 5nM+C','A 5nM+C','D 10nM+C','A 10nM+C', 'LF 1nM+C','HF 1nM+C','LF 5nM+C','HF 5nM+C','LF 10nM+C','HF 10nM+C'}

% i_donly = [25 27 29];
% i_aonly = [26 28 30 ];
i_donly = [3 5 7 ];
i_aonly = [4 6 8 ];

i_remove = [25 27 29 26 28 30 ];

dd = dd-mean(dd(i_bg));
aa = aa-mean(aa(i_bg));
da = da-mean(da(i_bg));

i_cells_only = [21 22 23 24];
i_cells = 20:36;

dd(i_cells) = dd(i_cells)-mean(dd(i_cells_only));
aa(i_cells) = aa(i_cells)-mean(aa(i_cells_only));
da(i_cells) = da(i_cells)-mean(da(i_cells_only));




leak = mean(da(i_donly)./dd(i_donly));
dir = mean(da(i_aonly)./aa(i_aonly));
da = da - leak.*dd - dir.*aa;

cur_fig = figure(2); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 20 15 ], 'PaperSize', [20 15] );

subplot(2, 1, 1)
plot(1:N_sample, dd, 'g.-', ...
    1:N_sample, aa, 'r.-', ...
    1:N_sample, da, 'b.-')
grid on
set(gca, 'XTick', 1:N_sample, 'XTickLabel', titles)
xtickangle(45)
ylabel('Median intensity')

subplot(2, 1, 2)
i_plot = setdiff(1:N_sample, [i_bg i_donly i_aonly i_cells_only i_remove])
plot(i_plot, da(i_plot)./(da(i_plot)+dd(i_plot)), 'k.')
grid on
set(gca, 'XTick', 1:N_sample, 'XTickLabel', titles)
xtickangle(45)
ylabel('FRET efficiency')

print(cur_fig, '-dpdf', [path_out filesep prefix_out '_analysis.pdf']); %save figure




%% -------------






%%
figure(3), clf
for i=1:N_sample

    
    tmp2 = tmp(:,:,i,1);
    tmp2 = tmp2-median(tmp2(:));
    
    
    tmp3 = tmp(:,:,i,1);
    tmp3 = tmp3-median(tmp3(:));
    
    
    tmp4 = tmp(:,:,i,1);
    tmp4 = tmp4-median(tmp4(:));
    
%     
     subplot(1, 3, 1)
     imagesc(tmp2), colorbar, axis image
     subplot(1, 3, 2)
     imagesc(tmp3), colorbar, axis image
     subplot(1, 3, 3)
     imagesc(tmp4), colorbar, axis image
%     
%     [p_tmp, ~] = calculate_ration_of_areas(tmp2, tmp4, 'display', 'on')
%     p(i) = p_tmp(1);
%     
    pause
end

%% save data
close all
disp('Saving data...')
save([path_out prefix_out '_data.mat'])
disp('data saved...')


%% write corrected images
disp('Writing images')

da_cor = (gelData_raw.images{3}-mean(da(i_bg))) ...
    - leak.* (gelData_raw.images{1}-mean(dd(i_bg))) - dir.*(gelData_raw.images{2}-mean(aa(i_bg))) + mean(da(i_bg));

t = Tiff([path_out filesep 'da_cor+bg.tif'],'w');
t.setTag('Photometric',Tiff.Photometric.MinIsWhite);
t.setTag('BitsPerSample',16);
t.setTag('SampleFormat',Tiff.SampleFormat.UInt);
t.setTag('ImageLength',size(da_cor,1));
t.setTag('ImageWidth',size(da_cor,2));
t.setTag('SamplesPerPixel',1);
t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
t.write( uint16(da_cor)  );
t.close();
 
disp('Done.')



