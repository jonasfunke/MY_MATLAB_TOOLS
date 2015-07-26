function [ img_out ] = my_imagetransform(img,  mirror, angle )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if mirror
        img_out = flip(imrotate(img, angle, 'crop'),1);
    else
        img_out = imrotate(img, angle, 'crop');
    end

end

