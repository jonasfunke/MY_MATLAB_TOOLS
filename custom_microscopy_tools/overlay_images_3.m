%% create composite RGB image
close all, clear all, clc

%path_in = '/Users/jonasfunke/Dropbox (DIETZ LAB)/PLECTONIC/LAB/Data/2020-09-16_evos_m7000_test/killing2.2020-09-17-16-10-29/'

path_in = '/Users/jonasfunke/Dropbox (DIETZ LAB)/PLECTONIC/LAB/Data/2020-09-16_evos_m7000_test/nalm_cb_all.2020-09-17-12-16-05/';
%%
stem = '*_A02f00';

files_fl = dir([path_in stem 'd2.TIF']);
files_bf = dir([path_in stem 'd4.TIF']);
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

