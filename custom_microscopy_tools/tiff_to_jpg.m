close all, clear all, clc

path_in = '/Users/jonasfunke/Dropbox (Personal)/Plectonic_Experiments/data_CellCounter/2020-03-13_killing_assay_28h/RGB_images/';

files_in = dir([path_in '*.tiff']);

path_out = [path_in];

%
format_out = 'jpg';
bit_depth = 8;
for i= 1:length(files_in) %[1:7 9:10] %
    img = imread([files_in(i).folder filesep files_in(i).name]);
    file_out = [files_in(i).folder filesep files_in(i).name(1:end-4) 'jpg'];
    imwrite(uint8( img ), file_out, format_out,'BitDepth',bit_depth);

end
disp('done')
%%