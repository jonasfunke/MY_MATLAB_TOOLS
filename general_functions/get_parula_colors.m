function [cc] = get_parula_colors(N_colors)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    tmp = parula();
    cc = tmp(1:floor(size(tmp,1)/(N_colors-1)):end,:);%varycolor(3);
    
end

