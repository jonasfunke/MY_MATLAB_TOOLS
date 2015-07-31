function [] = img_to_tif( )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    [fname, pname] = uigetfile('*.img', 'Select img file');
    
    path_out = [pname  fname(1:end-4) '_tif'];
    mkdir(path_out)
    
    img = ReadImagic([pname filesep fname]); % read imagic file
    for i=1:size(img,3)
        imwrite(uint16(img(:,:,i)), [path_out filesep fname(1:end-4) sprintf('_%04i.tif', i)]);
    end
    disp([fname ' exported to tif.'])

    
end

