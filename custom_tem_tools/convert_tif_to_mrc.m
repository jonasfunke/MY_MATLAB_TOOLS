function [  ] = convert_tif_to_mrc( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
   [pname] = uigetdir(cd, 'Select directory with TIF files.');
    
    files = dir([pname filesep '*.TIF']);
    path_out = [pname '_mrc'];
    mkdir(path_out)
    N_img = length(files);
    tmp = imread([pname filesep files(1).name]);
   
    %img = zeros(size(tmp,1), size(tmp,2), N_img, class(tmp));
    
    for i=1:N_img
        cur_img = imread([pname filesep files(i).name]);
        
        tmp = double(cur_img(1:2048, 1:2048));
        tmp = tmp-min(tmp(:));
        
        file_out = [path_out filesep files(i).name(1:end-4) '.mrc'];
        
        

        WriteMRC(single(tmp*(2^15-1)./max(tmp(:))), 1, file_out)
        
        
    end
   
   
    disp(['Done'])

end