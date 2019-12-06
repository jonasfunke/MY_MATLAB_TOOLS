%% 
close all, clear all, clc

%path_in = '/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_FACS/2019-12-02_PositioningScreen/';

%path_in = '/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_FACS/2019-12-03-Facs/Positioning/';
%path_in = '/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_FACS/2019-12-04_LGv7_ConfigurationScreen/'

path_in = '/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_FACS/2019-08-28_PositioningScreen/Data/'

files = dir([path_in '*.fcs']);

%%

for i=1:length(files)-5
    
  
    parts = split(files(i).name(1:end-4), ["-", "+"]);
    index = parts{end};
    
%     if strcmp(index, "0")
%         name = 'CBonly';
%     else
%         name = index;
%     end
    name = index;
    new_name = [ name '-' files(i).name(12:end)];
    file_in = [files(i).folder filesep files(i).name];
    mkdir([files(i).folder filesep 'renamed'])
    file_out = [files(i).folder filesep 'renamed' filesep new_name];
    copyfile(file_in, file_out)
    
end

%%