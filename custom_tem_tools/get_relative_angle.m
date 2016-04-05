function [ dalpha ] = get_relative_angle(varargin)
%Load data from imageJ(aquired with the arrow tool) and combutes relative
%angle

    if isempty(varargin)
        r_slice = 8;
        r_angle = 7;
    else
        r_angle = varargin{1};
        r_slice = varargin{2};
    end
    [fname, pname] =uigetfile('*.*', 'Select output file of imageJ.');
    data = dlmread([pname fname], '\t', 1); % load data

    % convert absolute angles to relative angles

    dalpha = zeros(size(data,1)/2,2);
    for i=1:2:size(data,1)-1
        if data(i,r_slice) == data(i+1,r_slice)
            dalpha((i+1)/2,1) = abs(data(i+1,r_angle)-data(i,r_angle));
            dalpha((i+1)/2,2) = data(i,r_slice); 
            if dalpha((i+1)/2,1) > 180
                dalpha((i+1)/2,1) = 360-dalpha((i+1)/2,1);
            end
        else
            disp(['Warning: Out of slice sync at ' num2str(i) ])
        end
    end


    if sum(diff(dalpha(:,2)==0)) > 0
        disp('Warning: some particles in wrong order')
        dalpha = sortrows(dalpha,2);
    end

    dlmwrite([pname fname(1:end-4) '_angles.txt'], dalpha, '\t') % write output

    cf = figure();
    xhist = 0:5:max(dalpha(:,1));
    n = hist(dalpha(:,1), xhist);
    bar(xhist, n)
    pause
    close(cf)
    pause(0.1)

    export_mat_file = questdlg('Would you like to combine the angles with images?', ...
        'Export to mat', ...
        'Yes','No','No');

    if strcmp(export_mat_file, 'Yes')
        [fname_img, pname_img] = uigetfile('*.img', 'Select corresponding img file');
        img = ReadImagic([pname_img filesep fname_img]); % read imagic file
        particles = img(:,:,dalpha(:,2));
        name = inputdlg({'Name of stack:'}, 'Name', 1, {fname});
        angles = dalpha(:,1);
        history = [{['Created on ' datestr(now, 'yyyy-mm-dd:')]}; ...
            {['from ' pname_img fname_img]}];
        save([pname fname(1:end-4) '_angles.mat'], 'angles', 'history', 'name', 'particles')

    end


end

