%% startup
clear all, close all, clc
dalpha = 5;

%% load images
%pname=uigetdir(data_directory,'Choose a folder with particle images.'); % get pathname
%tmp = dir([pname filesep '*.tif']);
%fnames = {tmp.name}; % list of filenames

pname = '/Users/jonasfunke/Documents/Figure_1/Class_averages_selection/cropped_200x200/img-files';
tmp = dir([pname filesep 'class_averages*.tif']);

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

%% align images
dalpha = 2;

alpha = 0:dalpha:359;
n_rot = length(alpha);
img_size = 200;
best_fit = zeros(img_size, img_size, n_img);

best_fit(:,:,1) = images(:,:,1);

for i=2:n_img
    % compute xcorrelation
    xcor_img = zeros(img_size, img_size, n_rot); 
    for r=1:n_rot % loop through rotations
        tmp = imrotate(images(:,:,i), alpha(r),  'bicubic','crop'); % rotate image
        tmp_xcor = normxcorr2(tmp, images(:,:,1)); % x-correlate
       % if i==2
       %     tmp_xcor = normxcorr2(tmp, images(:,:,i-1)); % x-correlate
       % else
       %     tmp_xcor = normxcorr2(tmp, best_fit(:,:,i-1)); % x-correlate
       % end
        xcor_img(:,:,r) = tmp_xcor(100:299, 100:299);%tmp_img(box_size/2+1:end-box_size/2, box_size/2+1:end-box_size/2);    
       % cor(r) = corr2(tmp, images(:,:,i-1));
    end

    % find maximum
    [max_cor, i_max_cor] = max(xcor_img(:))
    [k, l, r] = ind2sub(size(xcor_img), i_max_cor);
    
 %   [max_cor, r] = max(cor);

     best_fit(:,:,i) = imrotate(images(:,:,i), alpha(r), 'bicubic', 'crop');
   % rgb = zeros(200,200,3);
   % rgb(:,:,1) = images(:,:,i-1);
   % rgb(:,:,2) = best_fit;
   % imshow(uint16(rgb))
    
  %  pause
end


%%

img_mon = zeros(200,200,1,n_img);

for i=1:n_img
   % subplot(1,2,1)
   % imagesc(images(:,:,i-1)), axis image, colormap gray
    
    %subplot(1,2,2)
    img_mon(:,:,1,i) = images(:,:,i);
    
end

%%
close all
montage(uint16(img_mon(:,:,:,1:5)), 'Size', [1 5])

%% translate

alpha = 0; %:dalpha:359;
n_rot = length(alpha);
img_size = 200;

dxdy = zeros(n_img, 2);
for i=1:n_img
    % compute xcorrelation
    tmp_xcor = normxcorr2(best_fit(:,:,i), best_fit(:,:,1)); % x-correlate

    xcor_img = tmp_xcor(100:299, 100:299);%tmp_img(box_size/2+1:end-box_size/2, box_size/2+1:end-box_size/2);    

    % find maximum
    [max_cor, i_max_cor] = max(xcor_img(:));
    [k, l] = ind2sub(size(xcor_img), i_max_cor);
    dxdy(i,:) = [ k-101, l-101];
    
end
%%

img_out = zeros(170,170,n_img);

for i=1:n_img
    J = imtranslate(best_fit(:,:,i), [dxdy(i,2) dxdy(i,1)],'FillValues',139);
    img_out(:,:,i) = J(15:184, 15:184);
    
    subplot(1,2,1)
    imagesc(best_fit(:,:,i)), colormap gray, axis image
    subplot(1,2,2)
    imagesc(J), colormap gray, axis image
    
    imwrite(uint16(img_out(:,:,i)), [pname filesep 'rotated_' num2str(i) '.tif'])
   % pause
end










