%% welcome message
display('----------------------------------------------------------------')
display('---------------- Do. Or do not. There is no try. ---------------')
display('----------------------------------------------------------------')

%% set data directory
data_dir = data_directory(); % retrive where most data is stored
matlab_dir = userpath;
matlab_dir = matlab_dir(1:end-1); % home matlab folder

%% set default plot styles
set(0, 'defaulttextinterpreter', 'none')
%set(0,'defaultAxesFontName', 'Times-Roman')
%set(0,'defaultTextFontName', 'Times-Roman')
set(0,'defaultlinelinewidth',2)
set(0,'DefaultLineMarkerSize',15)

%% include the following 
% add directories, the subdirecotries will be added
subdirs_to_add = {...
    fullfile(matlab_dir, 'MATLAB_TOOLBOX'), ...
    fullfile(matlab_dir, 'MY_MATLAB_TOOLS'), ...
    %fullfile(matlab_dir, 'Add yur directory here'), ...
    };
dirs_to_add = {...
    %fullfile(matlab_dir, 'Add yur directory here'), ...
    };

for i=1:length(subdirs_to_add)
    if exist(subdirs_to_add{i}) == 7 %is a folder 
        add_subdirectories(subdirs_to_add{i})
    end
end

for i=1:length(dirs_to_add)
    if exist(subdirs_to_add{i}) == 7 %is a folder 
        addpath(dirs_to_add{i})
        disp(['Added: ' dirs_to_add{i}])
    end
end

%% change to matlab directory
cd(matlab_dir)

%% set default figureposition to left screen
%if size(get(0,'MonitorPositions'),1)>1
%    p = [1920+1280-560/2, 420/2,  560, 420];
%    set(0, 'DefaultFigurePosition', p)
%end
%% clear all variable for a clean startup
clear

