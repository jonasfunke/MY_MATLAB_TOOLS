function add_subdirectories(directory, excluded )
%This function adds all folders in 'directory' to the searchpath
%   It exludes, folder that start with +, @, ., private or givenin excluded

     if ~exist('excluded', 'var')
         excluded = {};
     end
    % determine list of directories
    all_files = dir(directory);
    directories = {all_files([all_files.isdir]).name};
    excluded_start = {'.', '+', '@'};
    excluded_directories = [{'private'}, excluded];

    dirs2add = cell(0,1);
    for i=1:length(directories)
        if ~strcmp(directories{i}, excluded_directories) % dirs that are excluded 
            if ~strncmp(directories{i}, excluded_start, 1) % dirs that staart with
                dirs2add = [dirs2add, fullfile(directory, directories{i})];
            end
        end
    end
    dirs2add = [directory, dirs2add];
    % add directories
    for i=1:length(dirs2add)
        addpath(dirs2add{i})
        disp(['Added: ' dirs2add{i}])
    end


end

