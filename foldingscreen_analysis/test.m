

%%
close all, clear all, clc

parsed_data = parse_gel_info_simple('/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/JF_Plate_v3/gel_info_simple.txt');

check_parsed_gel_info(parsed_data);

%%
close all, clear all, clc

parsed_data = parse_gel_info('/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/JF_FSv7/gel_config_initial_folding_screen.txt');

check_parsed_gel_info(parsed_data);


%%
parsed_data = parse_gel_info_simple('/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/Fako_csv2-45deg-shortoligos/Fako_csv2-45deg-shortoligos_initialfoldingscreen_gel_info.txt');

%%


root_path = '/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/';

folders = dir(root_path);

discard = {'.', '..', 'TEMPLATE'};

i_discard = [];
for i=1:length(folders)
    if any(strcmp(discard, folders(i).name))
        disp(['Discarding folder ' folders(i).name])
        i_discard = [i_discard, i];
    end
end
folders(i_discard) = [];

for i=1:length(folders)
    if folders(i).isdir
        disp('---')
        txt_files = dir([folders(i).folder filesep folders(i).name filesep '*.txt']);
        tif_files = dir([folders(i).folder filesep folders(i).name filesep '*.tif']);
        
        if length(txt_files)==1
            [tmp, warnings] = parse_gel_info_simple([txt_files(1).folder filesep txt_files(1).name]);
                check_parsed_gel_info(tmp);
        else
            disp(['More or less than one txt file found in ' folders(i).name ' (' num2str(length(txt_files)) ' found).' ])
        end
        
        if length(tif_files)==1
            %[tmp, warnings] = parse_gel_info_simple([txt_files(1).folder filesep txt_files(1).name]);
            %    check_parsed_gel_info(tmp);
        else
            disp(['More or less than one tif file found in ' folders(i).name ' (' num2str(length(tif_files)) ' found).' ])
        end
       pause 
       pause 
       
    end
end

%% Test compute_profiles

root_path = '/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS';
dir_name = 'JF_Plate_v3';
gel_info_name = 'gel_info_simple.txt';
image_name = '2018-12-15_Platev3_initial-folding-screen_2_25um-[EtBr].tif';

compute_profiles(root_path, dir_name, gel_info_name, image_name)

%%
close all, clear all, clc
root_path = '/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS';

[path_selected] = uigetdir(root_path, 'Select an Directory with intial foling screen');
[~,dir_name,] = fileparts(path_selected);

txt_files = dir([path_selected filesep '*.txt']);
tif_files = dir([path_selected filesep '*.tif']);

if length(tif_files)==1 && length(txt_files)==1
    compute_profiles(root_path, dir_name, txt_files(1).name, tif_files(1).name)
else
    disp(['Error. Found ' num2str(length(tif_files)) ' tif files and ' num2str(length(txt_files)) ' txt files in ' dir_name])
end
%




