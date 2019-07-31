function [ I, areas ] = get_area_intensities(img, n_areas, varargin)
% Opens an GUI to select and integrate multiple areas in an image (in cell of images/channels)
%   Input: 
%           img = cell of images (matrices) (all must have the same size), required
%           n_areas = Number of areas to integrate, required
%           option same_size: determines whether all areas will have the
%           same size (as the first area)
%           option plot_factor: weight for plotting the different channels
%           Output: 
%               I = vector/matrix of mean intensities for each area
%               areas = positions of the selected areas
%           examples:
%               integrate_areas(images, 3) integrate 3 areas, not resizable
%               integrate_areas(images, 3, 'resizable', true) integrate 3 resizable areas
%               integrate_areas(images, 3, 'plot_factor', [0 0 1]) integrate 3 areas, only show the third channel              
% 


    % parse input variables
    p = inputParser;
    default_resizable = false;
    default_plot_factor = ones(1, max(size(img)));
    default_message = 'Select area';
    
    addRequired(p,'img',@iscell);
    addRequired(p,'n_areas',@isnumeric);

    addParameter(p,'resizable',default_resizable, @islogical);
    addParameter(p,'plot_factor',default_plot_factor, @(x) length(img)==length(x) ); % check that plot factor has correct size
    addParameter(p,'message',default_message, @isstr);

    parse(p, img, n_areas, varargin{:});
    resizable = p.Results.resizable;
    plot_factor = p.Results.plot_factor;
    message = p.Results.message;
    
    % init variables
    areas = zeros(n_areas, 4);
    %I = zeros(n_areas, max(size(img)));
    N_img = max(size(img));
    
    % calculate image  thah will be used to select areas
    img_show = zeros(size(img{1},1), size(img{1},2));
    for i=1:N_img
        img_show = img_show + plot_factor(i) .* img{i};
    end
    
    
    % select areas by hand
    cur_fig = plot_image_ui(img_show); % show image
    
    % Create back-button
    Back_button = uicontrol('Style', 'pushbutton', 'String', 'Back',...
        'Position', [600 20 50 20],... %location, values based on plot_image_gui
        'Callback', @back);     
    
    title({message, ['Resizable = ' num2str((resizable)) ', plot_weights = ' num2str(plot_factor)]})
    i = 1;
    my_rectangles = cell(n_areas,1); % cell that stores all rectangle objects
    index_tex = cell(n_areas,1); % cell that stores all text objects
    while i <= n_areas
        if i==1
            h = imrect; % new area
        else
            if ~resizable
                h = imrect(gca, double(pos)); % new unresizeable area
                setResizable(h, 0)
            else
                h = imrect; % new area
            end
        end
        wait(h);
        pos = int32(getPosition(h)); % get position of current area[xmin ymin width height]
        delete(h)
        

        areas(i,:) = pos; % set output

        if i==1
            I = zeros(pos(4)+1, pos(3)+1, n_areas, max(size(img)));
        end
        
        %integrate channels and devide by area
        for j= 1:max(size(img))
            I(:,:,i,j) = img{j}(pos(2):pos(2)+pos(4)  , pos(1):pos(1)+pos(3)); %integrate
        end

        % plot selected area
        my_rectangles{i} = rectangle('Position', areas(i,:), 'EdgeColor', 'r'); %plot all integrated areas
        index_tex{i} = text(areas(i, 1)+areas(i,3)/2, areas(i, 2)+areas(i,4)/2, num2str(i), ...
            'Color', 'r', 'VerticalAlignment', 'Middle', 'HorizontalAlignment', 'Center') ; % plot index of band
        i = i+1; % increment counter
        
    end
    close(cur_fig)
    
    % callback function for back-button
    function back(source,callbackdata)
        if i > 1
            my_rectangles{i-1}.delete; % delete previous rectangle
            index_tex{i-1}.delete; % delete number
            i = i-1; %decrement counter
        else
            disp('Can not go back. Already at first areas.')
        end
    end
    
    
end
