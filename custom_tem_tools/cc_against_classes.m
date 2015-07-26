function [ max_cor ] = cc_against_classes( image, references )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
           
    max_cor = zeros(size(references, 3), 1);
    for i=1:size(references, 3)
        tmp = normxcorr2(references(:,:,i), image);
        index = size(image,1)/2:size(image,1)+size(image,1)/2; % 100:300 for a 200x200 image
        max_cor(i) = max( max( tmp(index,index) ) ); % x-correlate
    end

end

