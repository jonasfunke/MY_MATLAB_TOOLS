%% startup
close all, clear all, clc

%% load gel data
gelData_raw = load_gel_image('data_dir', data_directory);

%% background correct data
gelData = background_correct_gel_image(gelData_raw, 'numberOfAreas', 4);

%% check for saturation
gelData = check_gel_saturation(gelData);

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
leak_dir  = calculate_corrections(gelData.images{1}, gelData.images{3}, gelData.images{2}, [path_out filesep prefix_out '_correction.txt']);
da_cor = gelData.images{3} - leak_dir(1,1).*gelData.images{1} - leak_dir(2,1).*gelData.images{2};%-leak_dir(1,2)-leak_dir(2,2); 

%% integrate bands
bandData = get_band_intensities(gelData);

%%
gamma =0.1;
gamma_calc =  bandData.intensities(end,3).*(1./0.5 - 1) ./  bandData.intensities(end,1) 
E = bandData.intensities(:,3) ./ (gamma_calc.*bandData.intensities(:,1) + bandData.intensities(:,3));
n_bands = size(bandData.intensities,1);
close all
plot(1:n_bands, E, 'b.-')


%%
scrsz = get(0,'screensize');
fig_dim =[15 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

subplot(1,2,1)
plot(1:n_bands, E, 'b.-', 1:n_bands, bla.E, 'r.-')
legend({'no-stack: 2001', 'stack: 0021'})
set(gca, 'YLim', [0.2 0.32], 'XLim', [0 20])
xlabel('Lane')
ylabel('Raw FRET efficiency')

subplot(1,2,2)
plot(c, E(1:end-1),'b.-', c,  bla.E(1:end-1), 'r.-')
legend({'no-stack: 2001', 'stack: 0021'})
set(gca, 'YLim',  [0.2 0.32], 'XLim', [0 2])
xlabel('Additional NaCl [M]')
ylabel('Raw FRET efficiency')

print(cur_fig, '-dtiff', '-r 500' , [path_out filesep prefix_out '_raw-fret.tif']); %save figure

%%
close all

y1 = mean(E(1:5));
x1 = mean(bla.E(1:5));
y2 = E(end);
x2 = bla.E(end);
a = (y2-y1)./ (x2-x1)
b = (x1*y2-x2*y1) / (x1-x2)
plot(1:n_bands, E, 'b.-', 1:n_bands, bla.E.*a+b, 'r.-')



%%
E1 = E-mean(E(16:17));
E1 = E1/ median(E1(1:5));
E2 = bla.E-mean(bla.E(16:17));
E2 = E2/ median(E2(1:5));


%%
cur_fig = figure();
plot(1:n_bands, E1, 'b.-', 1:n_bands, E2, 'r.-')
legend({'no-stack: 2001', 'stack: 0021'})
set(gca, 'YLim', [-0.1 1.2], 'XLim', [0 20])
xlabel('Lane')
ylabel('Normalized FRET efficiency')

print(cur_fig, '-dtiff', '-r 500' , [path_out filesep prefix_out '_fret.tif']); %save figure

%%
cur_fig = figure();
plot(c, E1(1:end-1), 'b.-', c, E2(1:end-1), 'r.-')
legend({'no-stack: 2001', 'stack: 0021'})
set(gca, 'YLim', [-0.1 1.2], 'XLim', [0 2])
xlabel('Additional NaCl [M]')
ylabel('Normalized FRET efficiency')

print(cur_fig, '-dtiff', '-r 500' , [path_out filesep prefix_out '_fret_conc.tif']); %save figure


%% save data
save([path_out prefix_out '_data.mat'])


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



%%