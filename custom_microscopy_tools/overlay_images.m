%% create composite RGB image
close all, clear all, clc

%path_in = '/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_CellCounter/2020-03-03_Recruiting_CellTrace_20h/';
%path_in = '/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_CellCounter/2020-03-02_Recruiting_CellTrace_1.5h/';
%path_in = '/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_CellCounter/2020-03-02_Recruiting_CellTrace_4h/';
path_in = '/Users/jonasfunke/Dropbox (Personal)/Plectonic_Experiments/data_CellCounter/2020-03-04_Recruiting_42h/';
files_BF = dir([path_in '*_BF.tiff']);
% 
% path_in ='/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_CellCounter/2020-03-02_Recruiting_CellTrace_0h/';
% files_BF = dir([path_in '*_BF.PNG']);

path_out = [path_in 'RGB_images'];
mkdir(path_out)
%%

% high pass filter
r_filter = 100; % pixel
f_filter = fspecial('gaussian', 2*r_filter , r_filter); % gaussian filter
 

for i= 1:length(files_BF) %[1:7 9:10] %
    disp([num2str(i) ' of ' num2str(length(files_BF))]);
    BF_in = rgb2gray(imread([files_BF(i).folder filesep files_BF(i).name]));
    R_in = rgb2gray(imread([files_BF(i).folder filesep files_BF(i).name(1:end-7) 'CY5.tiff']));
    G_in = rgb2gray(imread([files_BF(i).folder filesep files_BF(i).name(1:end-7) 'RFP.tiff']));
    
%     BF_in = imread([files_BF(i).folder filesep files_BF(i).name]);
%     R_in = rgb2gray(imread([files_BF(i).folder filesep files_BF(i).name(1:end-6) 'CY5.PNG']));
%     G_in = rgb2gray(imread([files_BF(i).folder filesep files_BF(i).name(1:end-6) 'RFP.PNG']));


    %
    tmp =  double(BF_in)-double(imfilter(BF_in, f_filter, 'same', 'replicate'));
    BF = (tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)));

    tmp =  double(R_in)-double(imfilter(R_in, f_filter, 'same', 'replicate'));
    R = (tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)));
    R = min(R/0.6,1);

    tmp =  double(G_in)-double(imfilter(G_in, f_filter, 'same', 'replicate'));
    G = (tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)));
    % rescalt
    G = min(G/0.6,1);


    RGB = zeros(size(BF,1), size(BF,2),3, 'uint8');
    %RGB(:,:,1) = uint8(R*(2^8-1));
    %RGB(:,:,2) = uint8(G*(2^8-1));
    %RGB(:,:,3) = uint8(BF*(2^8-1));

    RGB(:,:,1) = uint8((BF/2+R/2)*(2^8-1));
    RGB(:,:,2) = uint8((BF/2+G/2)*(2^8-1));
    RGB(:,:,3) = uint8(BF/2*(2^8-1));
    
    %alpha = min(1-(BF-mean(BF(:))+std(BF(:)))./(1-median(BF(:))+std(BF(:))), 1);
    %imwrite(RGB, [path_out filesep files_BF(i).name(1:end-7) '_RGB.png'], 'Alpha', alpha);
    imwrite(RGB, [path_out filesep files_BF(i).name(1:end-7) '_RGB.tiff']);
    
    
    
%     % find peaks
%     tmp=FastPeakFind(R);
%     xy_R = [tmp(1:2:end) tmp(2:2:end)];
%     
%     tmp=FastPeakFind(G);
%     xy_G = [tmp(1:2:end) tmp(2:2:end)];
%     
%     tmp=FastPeakFind(BF);
%     xy_BF = [tmp(1:2:end) tmp(2:2:end)];
%     
%     disp([num2str(size(xy_BF,1)) ' - ' num2str(size(xy_R,1)) ' - ' num2str(size(xy_G,1)) ' - ' num2str(size(xy_R,1)+size(xy_G,1))])
% 
%     
%     
%     d_min = 40;
%     d = zeros(size(xy_R,1), size(xy_G,1));
%     pairs = [];
%     for k=1:size(xy_R,1)
%         for l=1:size(xy_G,1)
%             d(k,l)=sqrt((xy_R(k,1)-xy_G(l,1)).^2 + (xy_R(k,2)-xy_G(l,2)).^2);
%             if d(k,l) < d_min
%                pairs = [pairs; [xy_R(k,1) xy_R(k,2) xy_G(l,1) xy_G(l,2)]];
%             end
%         end
%     end
%     figure(1), clf
%     subplot(2,2,1)
%     imagesc(BF), colormap gray, axis image, hold on
%     plot(xy_BF(:,1),xy_BF(:,2),'bo')
%     subplot(2,2,2)
%     imagesc(R), colormap gray, axis image, hold on
%     plot(xy_R(:,1),xy_R(:,2),'ro')
%     subplot(2,2,3)
%     imagesc(G), colormap gray, axis image, hold on
%     plot(xy_G(:,1),xy_G(:,2),'go')
%     %
%     %
%     %
%     %pause
%     %plot(pairs(:,1),pairs(:,2),'r+')
%     %plot(pairs(:,3),pairs(:,4),'gx')
%     
    
end

%%
tmp = imregionalmax(R);
%%

figure(1), clf
imagesc(R), axis image, colorbar

%%
p=FastPeakFind(BF);
figure(1), clf
imagesc(BF); hold on, colorbar
plot(p(1:2:end),p(2:2:end),'r+')

%%



figure(1), clf
subplot(1,2,1)
imagesc(BF), axis image, colorbar

subplot(1,2,2)
imagesc(tmp),  axis image, colorbar

%%
subplot(1,3,3)
imagesc(rgb2gray(G)),  axis image, colorbar