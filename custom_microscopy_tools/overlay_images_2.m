%% create composite RGB image
close all, clear all, clc

%path_in = '/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_CellCounter/2020-03-03_Recruiting_CellTrace_20h/';
%path_in = '/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_CellCounter/2020-03-02_Recruiting_CellTrace_1.5h/';
%path_in = '/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_CellCounter/2020-03-02_Recruiting_CellTrace_4h/';
%path_in = '/Users/jonasfunke/Dropbox (Personal)/Plectonic_Experiments/data_CellCounter/2020-03-04_Recruiting_42h/';
%path_in = '/Users/jonasfunke/Dropbox (Personal)/Plectonic_Experiments/data_CellCounter/2020-03-12_killing-assay_0h/';
%path_in = '/Users/jonasfunke/Dropbox (Personal)/Plectonic_Experiments/data_CellCounter/2020-03-12_killing_assay_6h/'
%path_in = '/Users/jonasfunke/Dropbox (Personal)/Plectonic_Experiments/data_CellCounter/2020-03-12_killing_assay_3h/'
path_in = '/Users/jonasfunke/Dropbox (Personal)/Plectonic_Experiments/data_CellCounter/2020-03-13_killing_assay_28h/'
files_BF = dir([path_in '*_BF.tiff'])
% 
% path_in ='/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_CellCounter/2020-03-02_Recruiting_CellTrace_0h/';
% files_BF = dir([path_in '*_BF.PNG']);

path_out = [path_in 'RGB_images_2'];
mkdir(path_out)
%

% high pass filter
r_filter = 100; % pixel
f_filter = fspecial('gaussian', 2*r_filter , r_filter); % gaussian filter
 
hist_g = [];
for i= 1:length(files_BF) %[1:7 9:10] %
    disp([num2str(i) ' of ' num2str(length(files_BF))]);
    BF_in = rgb2gray(imread([files_BF(i).folder filesep files_BF(i).name]));
    %R_in = rgb2gray(imread([files_BF(i).folder filesep files_BF(i).name(1:end-7) 'CY5.tiff']));
    G_in = rgb2gray(imread([files_BF(i).folder filesep files_BF(i).name(1:end-7) 'RFP.tiff']));
    
%     BF_in = imread([files_BF(i).folder filesep files_BF(i).name]);
%     R_in = rgb2gray(imread([files_BF(i).folder filesep files_BF(i).name(1:end-6) 'CY5.PNG']));
%     G_in = rgb2gray(imread([files_BF(i).folder filesep files_BF(i).name(1:end-6) 'RFP.PNG']));


    %
    tmp =  double(BF_in)-double(imfilter(BF_in, f_filter, 'same', 'replicate'));
    BF = (tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)));

%     tmp =  double(R_in)-double(imfilter(R_in, f_filter, 'same', 'replicate'));
%     R = (tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)));
%     R = min(R/0.6,1);

    tmp =  double(G_in)-double(imfilter(G_in, f_filter, 'same', 'replicate'));
    %G = (tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)));
    new_max = 30;
    G = (tmp-mean(tmp(:))-3*std(tmp(:)))/(new_max-mean(tmp(:))-3*std(tmp(:)));
    % rescalt
    histogram(G(:), 20)
    pause
    G = max(G,0);
    G = min(G,1);
    

    RGB = zeros(size(BF,1), size(BF,2),3, 'uint8');
    %RGB(:,:,1) = uint8(R*(2^8-1));
    %RGB(:,:,2) = uint8(G*(2^8-1));
    %RGB(:,:,3) = uint8(BF*(2^8-1));

    RGB(:,:,1) = uint8((BF/4)*(2^8-1));
    RGB(:,:,2) = uint8((BF/4+G/2)*(2^8-1));
    RGB(:,:,3) = uint8(BF/4*(2^8-1));
    
    %alpha = min(1-(BF-mean(BF(:))+std(BF(:)))./(1-median(BF(:))+std(BF(:))), 1);
    %imwrite(RGB, [path_out filesep files_BF(i).name(1:end-7) '_RGB.png'], 'Alpha', alpha);
    imwrite(RGB, [path_out filesep files_BF(i).name(1:end-7) '_RGB.tiff']);
    imwrite(uint8( RGB ), [path_out filesep files_BF(i).name(1:end-7) '_RGB.jpg'], 'jpg' ,'BitDepth',8);

    
end
