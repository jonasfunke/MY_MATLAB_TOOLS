function [] = update_peaks(peaks, min_height,  fig_handle, markerspec)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %imagesc( img , 'Parent', fig_handle);
    %axis(fig_handle, 'image');
    xy = peaks(peaks(:,4)>=min_height, 1:2);

    
    children = get(fig_handle, 'children'); % get all previous plots on the cur_img
    delete(children(1:end-1)); % delete all previous plots (except for the last one = image)
    plot(xy(:,1), xy(:,2), markerspec, 'parent', fig_handle) % plot all the current peaks
end