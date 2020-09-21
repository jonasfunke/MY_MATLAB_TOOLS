%% create composite RGB image
close all, clear all, clc

path_in = '/Users/jonasfunke/Dropbox (DIETZ LAB)/PLECTONIC/LAB/Data/2020-09-16_evos_m7000_test/killing2.2020-09-17-16-10-29/'

%%
files_fl = dir([path_in '*_B02f00d0.TIF']);
files_bf = dir([path_in '*_B02f00d4.TIF']);
path_out = [files_bf(1).folder filesep files_bf(1).name(1:end-4) '/'];
mkdir(path_out)



figure(1); clf
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
    imwrite(imresize(RGB, 0.5), [path_out files_bf(i).name(1:end-4) '.png'], 'png' ,'BitDepth',8);

end
