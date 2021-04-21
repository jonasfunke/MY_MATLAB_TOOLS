%% create composite RGB image
close all, clear all, clc


path_in = '/Users/jonasfunke/Dropbox (DIETZ LAB)/PLECTONIC/LAB/Data/data_Microscope/SpotDetector.V4.1_03-15-21_06;54;24/';
%%
stem = '*';

files_fl = dir([path_in stem 'd1.TIFF']);
files_bf = dir([path_in stem 'd2.TIFF']);
%path_out = [files_fl(1).folder '_out' filesep files_fl(1).name(1:end-4) '/'];
path_out = [path_in 'zoom_ins' filesep]
mkdir(path_out)
%%
f_lim = [200 1000];
delta = 100; %px


for i= 1:1 %length(files_fl) %[1:7 9:10] %
    img_fl = double(imread([files_fl(i).folder filesep files_fl(i).name]));
    img_bf = double(imread([files_bf(i).folder filesep files_bf(i).name]));
    
    img_fl = (img_fl-min(img_fl(:)))/(max(img_fl(:))-min(img_fl(:)));
    img_bf = (img_bf-min(img_bf(:)))/(max(img_bf(:))-min(img_bf(:)));
    
    img_rgb = zeros(size(img_fl,1), size(img_fl,2), 3, 'uint8');
    img_rgb(:,:,1) =  (2^8-1)*(0.5*img_bf+1*img_fl);
    img_rgb(:,:,2) = (2^8-1)*0.5*img_bf ;
    img_rgb(:,:,3) = (2^8-1)*0.5*img_bf ;
    more_areas = 1;
    xy = [];
    while more_areas
        figure(1), clf
        imshow(img_rgb), hold on
        %rectangle('Position',[delta,delta,size(img_fl,2)-2*delta,size(img_fl,1)-2*delta])
        % select point of interest
        if ~isempty(xy)
            for k=1:size(xy, 1)
                plot(xy(:,1), xy(:,2), 'gx')
            end
        end
        [x_selected,y_selected] = ginput(1);
        x_selected = round(x_selected); 
        y_selected = round(y_selected); 
        
        xy = [xy; [x_selected,y_selected] ];

        % create zoom in
        cur_fig = figure(2); clf
        cur_ylim = [max(y_selected-delta,1) min(y_selected+delta,size(img_fl,1))];
        cur_xlim = [max(x_selected-delta,1) min(x_selected+delta,size(img_fl,2))];
        subplot(1, 2, 1)
        imagesc(img_bf(cur_ylim(1):cur_ylim(2), cur_xlim(1):cur_xlim(2)))
        colormap gray, axis image

        subplot(1, 2, 2)
        imagesc(img_fl(cur_ylim(1):cur_ylim(2), cur_xlim(1):cur_xlim(2)))
        colormap gray, axis image
        
        floc_out = [path_out files_fl(i).name(1:end-5) '_area' num2str(size(xy,1)) '.png'];
    
        set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
            'PaperPosition', [0 0 20 10 ], 'PaperSize', [20 10] );
        print(cur_fig, '-dpng', floc_out); %save figure

    
        answer = questdlg('Select more areas?','More areas',...
                  'Yes','No','No');
      if strcmp(answer, 'Yes')
          more_areas = 1;
      else
          more_areas = 0;
      end
    end
    
    
end


