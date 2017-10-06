function [  ] = crop_tem_micrographs( filter_bool, varargin )
%crop micrographs

    % parse input variables
    p = inputParser;
    default_filter = false;
    default_format = 'tif';
        
    addRequired(p,'filter_bool',@islogical);

    addParameter(p,'format',default_format, @isstr);
    
    parse(p, filter_bool,  varargin{:});
    format_out = p.Results.format;
    
  
    [ pname] =uigetdir('Select directory.');
    
    if filter_bool
        path_out = [pname '_cropped_filtered'];
    else
        path_out = [pname '_cropped'];
    end
    mkdir(path_out)
    
    files = dir([pname filesep '*.TIF']);
    
    % filter stuff
    r_filter = 50;
    f_filter = fspecial('gaussian', 2*r_filter , r_filter); % gaussian filter
   
    if strcmp('jpg', format_out)
        bit_depth = 8;
    else
        bit_depth = 16;
    end
       
    if filter_bool
        disp(['Filtering with high-pass gaussian radius of ' num2str(r_filter) ' px'])
    else
        disp('No filter')
    end
    disp(['Output: ' num2str(bit_depth) ' bit ' num2str(format_out)])
    
    for i=1:length(files)
        if strcmp('tif', format_out)
            file_out = [path_out filesep files(i).name];
        else
            file_out = [path_out filesep files(i).name(1:end-4) '.' format_out];
        end
        tmp = imread([pname filesep files(i).name]);
        if filter_bool
            tmp2 =  double(tmp(1:2048, 1:2048))-double(imfilter(tmp(1:2048, 1:2048), f_filter, 'same', 'replicate'));
            tmp2 = tmp2-min(tmp2(:));
            if bit_depth == 8
                imwrite(uint8( tmp2*(2^bit_depth-1)./max(tmp2(:)) ), file_out, format_out,'BitDepth',bit_depth);
            else
                imwrite(uint16( tmp2*(2^bit_depth-1)./max(tmp2(:)) ), file_out, format_out,'BitDepth',bit_depth);
            end
        else
            imwrite(tmp(1:2048, 1:2048), file_out, format_out,'BitDepth',bit_depth);
        end
    
    end
   
    disp('Micrographs cropped.')

end
