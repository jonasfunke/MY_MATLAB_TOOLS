%% create composite RGB image
close all, clear all, clc

path_in = '/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_CellCounter/2020-03-03_Recruiting_CellTrace_20h/';

files_BF = dir([path_in '*_BF.tiff']);

%%

% high pass filter
r_filter = 100; % pixel
f_filter = fspecial('gaussian', 2*r_filter , r_filter); % gaussian filter
 


for i=1:length(files_BF)
    i
    BF_in = rgb2gray(imread([files_BF(i).folder filesep files_BF(i).name]));
    R_in = rgb2gray(imread([files_BF(i).folder filesep files_BF(i).name(1:end-7) 'CY5.tiff']));
    G_in = rgb2gray(imread([files_BF(i).folder filesep files_BF(i).name(1:end-7) 'RFP.tiff']));


    %
    tmp =  double(BF_in)-double(imfilter(BF_in, f_filter, 'same', 'replicate'));
    BF = (tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)));

    tmp =  double(R_in)-double(imfilter(R_in, f_filter, 'same', 'replicate'));
    R = (tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)));

    tmp =  double(G_in)-double(imfilter(G_in, f_filter, 'same', 'replicate'));

    G = (tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)));



    RGB = zeros(size(BF,1), size(BF,2),3, 'uint8');
    RGB(:,:,1) = uint8(R*(2^8-1));
    RGB(:,:,2) = uint8(G*(2^8-1));
    %RGB(:,:,3) = uint8(BF*(2^8-1));

    alpha = min(1-(BF-mean(BF(:))+std(BF(:)))./(1-median(BF(:))+std(BF(:))), 1);
    imwrite(RGB, [files_BF(i).folder filesep files_BF(i).name(1:end-7) '_RGB.png'], 'Alpha', alpha);
    imwrite(RGB, [files_BF(i).folder filesep files_BF(i).name(1:end-7) '_RGB.tiff']);
end
%%
figure(1), clf
imagesc(RGB), axis image, colorbar

%%
%%



figure(1), clf
subplot(1,2,1)
imagesc(BF), axis image, colorbar

subplot(1,2,2)
imagesc(tmp),  axis image, colorbar

%%
subplot(1,3,3)
imagesc(rgb2gray(G)),  axis image, colorbar