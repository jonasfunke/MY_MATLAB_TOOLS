%% create composite RGB image
close all, clear all, clc


path_in = '/Users/jonasfunke/Dropbox (DIETZ LAB)/PLECTONIC/LAB/Data/data_Microscope/Jonas/2021-04-17_NALM-clustering-s8hb_8bit/';
%%
wells = {'H08', 'H09', 'H10', 'H11', 'H12'};
position = {'f01'};
time_points = cell(39,1);
for i=1:39
    time_points{i} = sprintf('%03i', i);
end

h = fspecial('average', 50);

for i=1:length(time_points) % loop through time points
   
    cur_fig =figure(1); clf
    
    for w=1:length(wells) % loop through wells
        
        stem = ['*'  time_points{i} wells{w} position{1}];
        files_fl = dir([path_in stem 'd1.TIF']);
        files_bf = dir([path_in stem 'd2.TIF']);
        
        img_bf = imresize(imread([files_bf(1).folder filesep files_bf(1).name]), 0.5);
        %img_fl = imread([files_fl(1).folder filesep files_fl(1).name]);
        
        cur_img = double(img_bf)-double(imfilter(img_bf, h, 'symmetric'));
        
        subplot(1,length(wells),w)
        imagesc(cur_img), axis image, colormap gray
        title([wells{w} ', t=' num2str((i-1)*0.5) 'h'])
        
        
    end
    set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
        'PaperPosition', [0 0 80 20 ], 'PaperSize', [80 20 ] );
    print(cur_fig, '-dpng',    [files_bf(1).folder filesep 'out_' time_points{i} '.png']); %save figure

end






%%
for w=1 %:length(wells) % loop through wells
    stem = ['*' wells{w} position{1}];
    files_fl = dir([path_in stem 'd1.TIF']);
    files_bf = dir([path_in stem 'd2.TIF']);
    
end

%%
files_fl = dir([path_in stem 'd1.TIF']);
files_bf = dir([path_in stem 'd2.TIF']);
path_out = [files_bf(1).folder '_out' filesep files_bf(1).name(1:end-4) '/'];
mkdir(path_out)



figure(1); clf
img_bf = imread([files_bf(1).folder filesep files_bf(1).name]);


v = VideoWriter( [path_out files_bf(1).name(1:end-6) ], 'MPEG-4');
v.FrameRate = 2;
open(v)
RGB = zeros(size(img_bf,1), size(img_bf,2),3, 'uint8');

for i= 1:length(files_bf) %[1:7 9:10] %
    disp([num2str(i) ' of ' num2str(length(files_bf))]);
    img_bf = imread([files_bf(i).folder filesep files_bf(i).name]);
    img_fl = imread([files_fl(i).folder filesep files_fl(i).name]);
    
    tmp = double(img_bf);
    bf_norm = (tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)));
    tmp = double(img_fl);
    fl_norm = (tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)));
    

    RGB(:,:,1) = uint8((bf_norm/1)*(2^8-1));
    RGB(:,:,2) = uint8((bf_norm/1 + (fl_norm-0.1)/1)*(2^8-1));
    RGB(:,:,3) = uint8((bf_norm/1)*(2^8-1));
 
    %imshow(RGB)
    %pause
    img_out = imresize(RGB, 0.5);
    imwrite(img_out, [path_out files_bf(i).name(1:end-6) '.png'], 'png' ,'BitDepth',8);
    writeVideo(v,img_out)
end
close(v)
disp('done')



%%


files_fl = dir([path_in  '*d2.TIF']);
files_bf = dir([path_in  '*d4.TIF']);

%%
path_out = [files_bf(1).folder '_out' filesep files_bf(1).name(1:end-4) '/'];
mkdir(path_out)



img_bf = imread([files_bf(1).folder filesep files_bf(1).name]);

RGB = zeros(size(img_bf,1), size(img_bf,2),3, 'uint8');
for i= 1:length(files_bf) %[1:7 9:10] %
    disp([num2str(i) ' of ' num2str(length(files_bf))]);
    img_bf = imread([files_bf(i).folder filesep files_bf(i).name]);
    img_fl = imread([files_fl(i).folder filesep files_fl(i).name]);
    
    tmp = double(img_bf);
    bf_norm = (tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)));
    tmp = double(img_fl);
    fl_norm = (tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)));
    

    RGB(:,:,1) = uint8((bf_norm/1)*(2^8-1));
    RGB(:,:,2) = uint8((bf_norm/1 + (fl_norm-0.1)/1)*(2^8-1));
    RGB(:,:,3) = uint8((bf_norm/1)*(2^8-1));
 
    %imshow(RGB)
    %pause
    img_out = RGB;
    imwrite(img_out, [path_out files_bf(i).name(1:end-6) '.png'], 'png' ,'BitDepth',8);
    imwrite(uint8(bf_norm*(2^8-1)), [path_out files_bf(i).name(1:end-6) '_BF.png'], 'png' ,'BitDepth',8);
    imwrite(uint8(fl_norm*(2^8-1)), [path_out files_bf(i).name(1:end-6) '_fl.png'], 'png' ,'BitDepth',8);
end
disp('done')

