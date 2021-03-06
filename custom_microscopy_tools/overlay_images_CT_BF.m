%% 
close all, clear all, clc



%path_in = '/Volumes/Jonas/2021-04-19_killing-A20-splenocytes/';
%path_in = '/Volumes/Jonas/2021-04-19_killing-A20-splenocytes_sytox/';
path_in = '/Volumes/Jonas/2021-04-22_killing-molm-s8hb-images/';

export_images = 1;
export_movie = 1;
framerate = 2;

lim = [0 255; ...
    0 100];

%stems = {'B02f01' 'B03f01' 'B04f01' 'B05f01' 'B06f01' 'B07f01' 'C02f01' 'C03f01' 'C04f01' 'C05f01' 'C06f01' 'C07f01'} ;
stems = {'E01f01'} ;
%

figure(1), clf
for k=1:length(stems)
    files_b = dir([path_in '*' stems{k} 'd1.TIF']);
    files_bf = dir([path_in '*' stems{k} 'd2.TIF']);
    %files_r = dir([path_in '*' stems{k} 'd3.TIF']);
    %files_g = dir([path_in '*' stems{k} 'd4.TIF']);


    path_out = [files_bf(1).folder '_' stems{k} filesep];
    mkdir(path_out)


    img_bf = imread([files_bf(1).folder filesep files_bf(1).name]);
    scale = 2;

    RGB = zeros(size(img_bf,1)/scale, size(img_bf,2)/scale, 3, 'uint8');
    img = zeros(size(img_bf,1)/scale, size(img_bf,2)/scale, 4, 'double');

    if export_movie   
        v = VideoWriter( [path_out files_bf(1).name(1:end-6) '_movie'], 'MPEG-4');
        v.FrameRate = framerate;
        open(v)
    end
    
    for i= 1:length(files_bf) %[1:7 9:10] %
        disp([num2str(i) ' of ' num2str(length(files_bf))]);
        img(:,:,1) = imresize(imread([files_bf(i).folder filesep files_bf(i).name]), 1/scale);
        %img(:,:,2) = imresize(imread([files_r(i).folder filesep files_r(i).name]), 1/scale);
        %img(:,:,3) = imresize(imread([files_g(i).folder filesep files_g(i).name]), 1/scale);
        img(:,:,2) = imresize(imread([files_b(i).folder filesep files_b(i).name]), 1/scale);
        
        tmp = img(:,:,2);
        histogram(tmp+1, logspace(0,3), 'DisplayStyle', 'stairs'), hold on
        set(gca, 'XLim', [1 256], 'xscale', 'log')
        %xline( prctile(tmp(:), 25), 'k')
        pause
        
        for j=1:2
            tmp = img(:,:,j);
            if j>1
                cur_min = prctile(tmp(:), 35);%  lim(j,1);
                cur_max = prctile(tmp(:), 98);%lim(j,2);
            else
                cur_min = lim(j,1);
                cur_max = lim(j,2);
            end
            
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
        RGB(:,:,1) = uint8((img(:,:,1)/2*(2^8-1) ));
        RGB(:,:,2) = uint8((img(:,:,1)/2*(2^8-1) ));
        RGB(:,:,3) = uint8( (img(:,:,1)/2 + img(:,:,2)/2)*(2^8-1));


        %imshow(RGB)
        %pause
        %img_out = RGB;
        if export_images
            %imwrite(img_out, [path_out files_bf(i).name(1:end-6) '.png'], 'png' ,'BitDepth',8);
            %imwrite(uint8(bf_norm*(2^8-1)), [path_out files_bf(i).name(1:end-6) '_BF.png'], 'png' ,'BitDepth',8);
            %imwrite(uint8(fl_norm*(2^8-1)), [path_out files_bf(i).name(1:end-6) '_fl.png'], 'png' ,'BitDepth',8);
            imwrite(RGB, [path_out files_bf(i).name(1:end-6) '.jpg'], 'jpg' ,'BitDepth',8);
        end
        if export_movie
            writeVideo(v,RGB)
        end
        
    end
    disp('done')
end
close(v)

