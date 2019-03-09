

%%
close all, clear all, clc

parsed_data = parse_gel_info_simple('/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/JF_Plate_v3/gel_info_simple.txt');

check_parsed_gel_info(parsed_data)


%%
parsed_data = parse_gel_info_simple('/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/Fako_csv2-45deg-shortoligos/Fako_csv2-45deg-shortoligos_initialfoldingscreen_gel_info.txt');

%%


root_path = '/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/'

folders = dir(root_path);

%%

for i=1:length(folders)
    if folders(i).isdir
        txt_files = dir([folders(i).folder filesep folders(i).name filesep '*.txt']);
        if isempty(txt_files)
            disp(['No txt files in ' folders(i).name])
        else
            if length(txt_files)>1
                disp(['More than one txt file in ' folders(i).name])
            else
                % parse
                disp('---')
                [tmp, warnings] = parse_gel_info_simple([txt_files(1).folder filesep txt_files(1).name]);
                check_parsed_gel_info(tmp);

            end
        end
        
    end
end

%% Test compute_profiles

root_path = '/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS';
dir_name = 'JF_Plate_v3';
gel_info_name = 'gel_info_simple.txt';
image_name = '2018-12-15_Platev3_initial-folding-screen_2_25um-[EtBr].tif';

compute_profiles(root_path, dir_name, gel_info_name, image_name)







