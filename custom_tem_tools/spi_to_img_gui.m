function [  ] = spi_to_img_gui(  )
% takes  a spi file and writes an img fle
%   Detailed explanation goes here
%% startup

path0 = cd;
data_location = '/Users/jonasfunke/Documents/FRET_STAGE/TEM images';

% get file
cd (data_location)
[fname, pname] = uigetfile('*.spi', 'Select a .spi file top be converted to img.');
cd(path0)

%% load spi file
images = readSPIDERfile([pname fname]);
disp(['Loaded: ' pname fname])

%% write img file
tmp = strsplit(pname, filesep);
prefix = tmp{end-2};
path_out = [pname prefix filesep];
mkdir(path_out);
file_out  = [path_out prefix '_ref_1'];
WriteImagic(images, file_out)
disp(['Writing ' file_out ' ... done'])
disp('Finished')


end

