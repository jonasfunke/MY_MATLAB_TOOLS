%% startup
close all, clear all, clc

%% load gel data
gelData_raw = load_gel_image('data_dir', data_directory, 'n_images', 3);

%% check for saturation
gelData_raw = check_gel_saturation(gelData_raw);

%% background correct data
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

%% shift images
%[dd_bg, da_x_min, da_y_min] = overlay_image(da_bg, dd_bg, 10);
%[aa_bg, aa_x_min, aa_y_min]= overlay_image(da_bg, aa_bg, 10);

%% leakage and direct-excitation correction factors
[leak_dir, leakdir_fig]  = calculate_leakage_directexcitation(gelData.images{1}, gelData.images{3}, gelData.images{2}, 'display', 'on');
print(leakdir_fig, '-dtiff', '-r 500' , [path_out filesep 'Leakage_DirEx_correction.tif']); %save figure

% write leakdir
save([path_out filesep 'Leakage_DirEx_data.txt'], 'leak_dir' ,'-ascii')

da_cor = gelData.images{3} - leak_dir(1,1).*gelData.images{1} - leak_dir(2,1).*gelData.images{2};%-leak_dir(1,2)-leak_dir(2,2); %correct image
gelData.images{4} = da_cor; % append to images

%% integrate bands
bandData = get_band_intensities(gelData);
pause(0.1)
close all
%% calculate intensities based on ratio
%{
I = zeros(18,3);
for i=1:18
    pos = bandData.positions(i,:);
    subDD = gelData.images{1}( pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3) );
    subAA = gelData.images{2}( pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3) );
    subDA = gelData.images{4}( pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3) );
    
    shiftDA = min(min(min(subDA))+0.000000001,0.000000001)
    shiftDD = min(min(min(subDD))+0.000000001,0.000000001)
        shiftAA = min(min(min(subAA))+0.000000001,0.000000001)

    E3(i) = sum(sum(( (subAA-shiftAA)./ sum(subAA(:)-shiftAA)).* (subDA-shiftDA) ./ ((subDA-shiftDA) + gamma_calc.*(subDD-shiftDD)))) ;
    
    bla =  (subDA-shiftDA) ./ ((subDA-shiftDA) + gamma_calc.*(subDD-shiftDD));
    prob = (subAA-shiftAA)./ sum(subAA(:)-shiftAA);
    [n, x] = histwc(bla(:), prob(:), 300);
    figure
    bar(x, n)
    title([num2str(i) ' weighted'])

    %figure
    %hist(bla(:), 100)
    
    %imagesc(   (subAA./ sum(subAA(:))).* (subDA+shiftDA) ./ ((subDA+shiftDA) + gamma_calc.*(subDD+shiftDD))) , colormap gray, axis image
    
    %pause
    p = calculate_ration_of_areas(subDD, subAA, 'display', 'off')
    I(i,1) = p(1);   
   % title(num2str(i))
    p = calculate_ration_of_areas(subDA, subAA, 'display', 'off')
    I(i,2) = p(1);
      %  title(num2str(i))

    p = calculate_ration_of_areas(subDD, subDA, 'display', 'off')
    I(i,3) = p(1);
      %  title(num2str(i))

end

%}


%%
DD_div_DA = zeros(size(bandData.intensities,1), 2);
DD_div_AA = zeros(size(bandData.intensities,1), 2);
DA_div_AA = zeros(size(bandData.intensities,1), 2);

for i=1:size(bandData.intensities,1)
    pos = bandData.positions(i,:);
    subDD = gelData.images{1}( pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3) );
    subAA = gelData.images{2}( pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3) );
    subDA = gelData.images{4}( pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3) );
    DD_div_DA(i,:) = calculate_ration_of_areas(subDD, subDA, 'display', 'on');
    pause
    close all
    DD_div_AA(i,:) = calculate_ration_of_areas(subDD, subAA, 'display', 'off');
    DA_div_AA(i,:) = calculate_ration_of_areas(subDA, subAA, 'display', 'off');
end

%%

i_gamma = 38;
E_soll = 0.5;
%gamma_calc =  bandData.intensities(i_gamma,4).*(1./0.5 - 1) ./  bandData.intensities(i_gamma,1) 
    
E_raw = 1./(1+DD_div_DA(:,1));
gamma_calc = (1-E_soll)./E_soll./DD_div_DA(i_gamma,1)
E = 1./(1+gamma_calc.*DD_div_DA(:,1));

gamma_calc_integrate =  bandData.intensities(i_gamma,4).*(1./E_soll - 1) ./  bandData.intensities(i_gamma,1) 
E_integrate = bandData.intensities(:,4) ./ (gamma_calc.*bandData.intensities(:,1) + bandData.intensities(:,4));


%% save data
close all
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
 
disp('Done.')

%% Plot 
n_bands = size(bandData.intensities,1);
cur_fig = figure;
subplot(3, 1, 1)
%plot(1:n_bands, bandData.intensities(:,2), 'r.-', 1:n_bands, bandData.intensities(:,1), 'g.-', 1:n_bands, bandData.intensities(:,4), 'b.-')
plot( 1:n_bands, gamma_calc.*bandData.intensities(:,1)./bandData.intensities(:,2), 'g.--', 1:n_bands, bandData.intensities(:,4)./bandData.intensities(:,2), 'b.--', ...
     1:n_bands, gamma_calc.*DD_div_AA(:,1), 'g.-', 1:n_bands, DA_div_AA(:,1), 'b.-')
xlabel('Lane'), ylabel('Normalized bandintensity')
legend({'gamma * D->D / A->A', 'D->A / A->A'}, 'location', 'best')
set(gca, 'XLim', [1 n_bands]);

subplot(3, 1, 2:3)
plot( 1:n_bands, E_integrate, 'k.--', 1:n_bands, E, 'k.-')
xlabel('Lane'), ylabel('FRET efficiency')

title({['gamma=' num2str(gamma_calc) ]})
legend({ 'FRET from intgration', 'FRET from scatterplot'})
set(gca, 'XLim', [1 n_bands]);

print(cur_fig, '-dtiff', '-r 500' , [path_out filesep 'FRET_normalized.tif']); %save figure
