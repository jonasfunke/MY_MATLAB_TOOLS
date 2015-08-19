function [  ] = convert_tif_to_img( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
   [pname] = uigetdir(cd, 'Select directory with TIF files.');
    
    files = dir([pname filesep '*.TIF']);

    N_img = length(files);
    tmp = imread([pname filesep files(1).name]);
   
    img = zeros(size(tmp,1), size(tmp,2), N_img, class(tmp));
    
    for i=1:N_img
        img(:,:,i) = imread([pname filesep files(i).name]);
    end
    
   
    fname = files(1).name(1:end-4);
    WriteImagic(img, [pname filesep fname])
    
   
    disp([fname ' converted to img.'])

end

