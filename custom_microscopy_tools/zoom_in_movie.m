%% 
close all, clear all, clc

path_in = '/Volumes/Jonas/2021-04-19_killing-A20-splenocytes/';


stems = { 'C06f01'} ;
delta = 300; %px
lim = [0 255; ...
    46 154; ...
    70 255; ...
   15 110];


for k=1:length(stems)
    files_b = dir([path_in '*' stems{k} 'd1.TIF']);
    files_bf = dir([path_in '*' stems{k} 'd2.TIF']);
    files_r = dir([path_in '*' stems{k} 'd3.TIF']);
    files_g = dir([path_in '*' stems{k} 'd4.TIF']);

    path_out = [files_bf(1).folder '_' stems{k} '_zoom5' filesep];
    mkdir(path_out)


    img_bf = imread([files_bf(1).folder filesep files_bf(1).name]);
    scale = 1;

    RGB = zeros(size(img_bf,1)/scale, size(img_bf,2)/scale, 3, 'uint8');
    img = zeros(size(img_bf,1)/scale, size(img_bf,2)/scale, 4, 'double');

    for i= 1:length(files_bf) %[1:7 9:10] %
        disp([num2str(i) ' of ' num2str(length(files_bf))]);
        img(:,:,1) = imresize(imread([files_bf(i).folder filesep files_bf(i).name]), 1/scale);
        img(:,:,2) = imresize(imread([files_r(i).folder filesep files_r(i).name]), 1/scale);
        img(:,:,3) = imresize(imread([files_g(i).folder filesep files_g(i).name]), 1/scale);
        img(:,:,4) = imresize(imread([files_b(i).folder filesep files_b(i).name]), 1/scale);

        if i==1
            % slect location
            figure(1), clf
            imagesc(img(:,:,1)), axis image, colormap gray, hold on
            [x_selected,y_selected] = ginput(1);
            x_selected = round(x_selected); 
            y_selected = round(y_selected); 
            cur_ylim = [max(y_selected-delta,1) min(y_selected+delta,size(img_bf,1))];
            cur_xlim = [max(x_selected-delta,1) min(x_selected+delta,size(img_bf,2))];
            img_zoom = zeros(cur_ylim(2)-cur_ylim(1)+1, cur_xlim(2)-cur_xlim(1)+1, 4, 'double');
            RGB = zeros(size(img_zoom,1), size(img_zoom,2), 3, 'uint8');
        end
        
        img_zoom = img(cur_ylim(1):cur_ylim(2), cur_xlim(1):cur_xlim(2), :);
        
        
        for j=1:4
            tmp = img_zoom(:,:,j);
            if j>1
                cur_min = prctile(tmp(:), 25);%  lim(j,1);
            else
                cur_min = lim(j,1);
            end
            cur_max = lim(j,2);
            %img(:,:,j) = (img(:,:,j)-min(min(img(:,:,j))))/(max(max(img(:,:,j)))-min(min(img(:,:,j))));
            img_zoom(:,:,j) = (img_zoom(:,:,j)-cur_min)/(cur_max-cur_min);
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
        RGB(:,:,1) = uint8((img_zoom(:,:,1)/2 + img_zoom(:,:,2)/2)*(2^8-1));
        RGB(:,:,2) = uint8((img_zoom(:,:,1)/2 + img_zoom(:,:,3)/2)*(2^8-1));
        RGB(:,:,3) = uint8((img_zoom(:,:,1)/2 + img_zoom(:,:,4)/2)*(2^8-1));


        %imshow(RGB)
        %pause
        %img_out = RGB;
        %imwrite(img_out, [path_out files_bf(i).name(1:end-6) '.png'], 'png' ,'BitDepth',8);
        %imwrite(uint8(bf_norm*(2^8-1)), [path_out files_bf(i).name(1:end-6) '_BF.png'], 'png' ,'BitDepth',8);
        %imwrite(uint8(fl_norm*(2^8-1)), [path_out files_bf(i).name(1:end-6) '_fl.png'], 'png' ,'BitDepth',8);
        imwrite(RGB, [path_out files_bf(i).name(1:end-6) '.png'], 'png' ,'BitDepth',8);
    end
    disp('done')
end

