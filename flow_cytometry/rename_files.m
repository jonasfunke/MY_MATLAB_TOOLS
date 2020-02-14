%% 
close all, clear all, clc

%path_in = '/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_FACS/2019-12-02_PositioningScreen/';

%path_in = '/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_FACS/2019-12-03-Facs/Positioning/';
%path_in = '/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_FACS/2019-12-04_LGv7_ConfigurationScreen/'

%path_in = '/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_FACS/2020-01-07-Concentration-Screen/NALM-CB-CD19/'

%path_in = '/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_FACS/2019-08-28_PositioningScreen/Data/'

path_in = '/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_FACS/2020-01-02-Valency-Screen-JURKAT/'
files = dir([path_in '*.fcs']);

%%
% make output dir
tmp = split(path_in, filesep);
dir_out = tmp{1};
for i=2:length(tmp)-2 % use / at the end
        dir_out = [dir_out filesep tmp{i}];
end
dir_out = [dir_out filesep [tmp{end-1} '_renamed']]
mkdir(dir_out)


%%

for i=1:length(files)-1
    
  
    [parts, matches] = split(files(i).name(1:end-4), ["-", "+"]);
    %index = parts{end};
    %time = parts{end-1};
    
%     if strcmp(index, "0")
%         name = 'CBonly';
%     else
%         name = index;
%     end
    new_name = parts{1};
    for j=2:length(parts)-2
        new_name = [new_name matches{j-1} parts{j}];
    end
    new_name = [new_name '-' parts{end} '-' parts{end-1} '.fcs'];
        
    %name = indirdex;
    %new_name = [ name '-' files(i).name(12:end)];
    file_in = [files(i).folder filesep files(i).name];

    file_out = [dir_out filesep new_name];
    copyfile(file_in, file_out)
    %disp(files(i).name)
    %disp(new_name)
    
end

%%