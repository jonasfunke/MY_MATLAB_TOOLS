function [] = normalize_img_stack( )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    [fname, pname] = uigetfile('*.img', 'Select img file');
    
    img = ReadImagic([pname filesep fname]); % read imagic file
    h = waitbar(0,'Images are being normalized...');
    img_out = zeros(size(img), 'uint16');
    for i=1:size(img,3)
        tmp = img(:,:,i)-min(img(:));
        img_out(:,:,i) =  (2^16-1).*tmp./max(tmp(:)) ;
        waitbar(i / size(img,3), h)
    end
    
    WriteImagic(img_out, [pname fname(1:end-4) '_normalized'])
    close(h)
    disp([fname ' normalized'])
    
end