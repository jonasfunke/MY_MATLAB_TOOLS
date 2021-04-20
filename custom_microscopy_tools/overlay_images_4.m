%% 
close all, clear all, clc



path_in = '/Volumes/Jonas/2021-04-19_killing-A20-splenocytes/';

stem = 'B02f01';

files_b = dir([path_in '*' stem 'd1.TIF']);
files_bf = dir([path_in '*' stem 'd2.TIF']);
files_r = dir([path_in '*' stem 'd3.TIF']);
files_g = dir([path_in '*' stem 'd4.TIF']);


path_out = [files_bf(1).folder '_' stem filesep];
mkdir(path_out)


%




lim = [0 255; ...
    46 154; ...
    70 255; ...
   15 110];

img_bf = imread([files_bf(1).folder filesep files_bf(1).name]);
scale = 2;

RGB = zeros(size(img_bf,1)/scale, size(img_bf,2)/scale, 3, 'uint8');
img = zeros(size(img_bf,1)/scale, size(img_bf,2)/scale, 4, 'double');

for i= 1:length(files_bf) %[1:7 9:10] %
    disp([num2str(i) ' of ' num2str(length(files_bf))]);
    img(:,:,1) = imresize(imread([files_bf(i).folder filesep files_bf(i).name]), 1/scale);
    img(:,:,2) = imresize(imread([files_r(i).folder filesep files_r(i).name]), 1/scale);
    img(:,:,3) = imresize(imread([files_g(i).folder filesep files_g(i).name]), 1/scale);
    img(:,:,4) = imresize(imread([files_b(i).folder filesep files_b(i).name]), 1/scale);
    
    for j=1:4
        tmp = img(:,:,j);
        if j>1
            cur_min = prctile(tmp(:), 25);%  lim(j,1);
        else
            cur_min = lim(j,1);
        end
        cur_max = lim(j,2);
        %img(:,:,j) = (img(:,:,j)-min(min(img(:,:,j))))/(max(max(img(:,:,j)))-min(min(img(:,:,j))));
        img(:,:,j) = (img(:,:,j)-cur_min)/(cur_max-cur_min);
    end
    
%     figure(1), clf
%     for j=1:4
%         subplot(2, 4, j)
%         tmp = img(:,:,j);
%         histogram(tmp(:))
%         
%         subplot(2, 4, j+4)
%         imagesc(img(:,:,j), [0 1]), axis image, colormap gray
%     end
    %pause
    RGB(:,:,1) = uint8((img(:,:,1)/2 + img(:,:,2)/2)*(2^8-1));
    RGB(:,:,2) = uint8((img(:,:,1)/2 + img(:,:,3)/2)*(2^8-1));
    RGB(:,:,3) = uint8((img(:,:,1)/2 + img(:,:,4)/2)*(2^8-1));
    
    
    %imshow(RGB)
    %pause
    %img_out = RGB;
    %imwrite(img_out, [path_out files_bf(i).name(1:end-6) '.png'], 'png' ,'BitDepth',8);
    %imwrite(uint8(bf_norm*(2^8-1)), [path_out files_bf(i).name(1:end-6) '_BF.png'], 'png' ,'BitDepth',8);
    %imwrite(uint8(fl_norm*(2^8-1)), [path_out files_bf(i).name(1:end-6) '_fl.png'], 'png' ,'BitDepth',8);
    imwrite(RGB, [path_out files_bf(i).name(1:end-6) '.jpg'], 'jpg' ,'BitDepth',8);
end
disp('done')

%%


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


path_out = [files_bf(1).folder '_out' filesep files_bf(1).name(1:end-4) '/'];
mkdir(path_out)


