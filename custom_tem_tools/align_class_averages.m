%% startup
clear all, close all, clc
dalpha = 5;

%% load images
pname=uigetdir(data_directory,'Choose a folder with particle images.'); % get pathname
tmp = dir([pname filesep '*.tif']);
fnames = {tmp.name}; % list of filenames
n_img = size(fnames,2);

% read images
tmp = imread([pname filesep fnames{1}]);
images = zeros(size(tmp,1), size(tmp,2), n_img);
for i=1:n_img
    images(:,:,i) = imread([pname filesep fnames{i}]);
end

% make output dir
path_out = [pname filesep datestr(now, 'yyyy-mm-dd_HH-MM') '_particles']; % output folder
mkdir(path_out)
%%
images(:,:,2) =  flipdim(images(:,:,2),1);

%% align images
dalpha = 5;

alpha = 0:dalpha:359;
n_rot = length(alpha);
img_size = 201;

for i=2:n_img
    % compute xcorrelation
    xcor_img = zeros(img_size, img_size, n_rot); 
    for r=1:n_rot % loop through rotations
        tmp = imrotate(images(:,:,i), alpha(r), 'crop'); % rotate image
        tmp_xcor = normxcorr2(tmp, images(:,:,i-1)); % x-correlate
        xcor_img(:,:,r) = tmp_xcor(100:300, 100:300);%tmp_img(box_size/2+1:end-box_size/2, box_size/2+1:end-box_size/2);    
        cor(r) = corr2(tmp, images(:,:,i-1));
    end

    % find maximum
    [max_cor, i_max_cor] = max(xcor_img(:))
    [k, l, r] = ind2sub(size(xcor_img), i_max_cor);
    
 %   [max_cor, r] = max(cor);
    best_fit = imrotate(images(:,:,i), alpha(r), 'crop');

    rgb = zeros(201,201,3);
    rgb(:,:,1) = images(:,:,i-1);
    rgb(:,:,2) = best_fit;
    imshow(uint16(rgb))
    
    

    pause
end
