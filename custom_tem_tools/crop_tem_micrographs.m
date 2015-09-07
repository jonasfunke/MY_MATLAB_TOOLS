function [  ] = crop_tem_micrographs( filter_bool )
%crop micrographs

    [ pname] =uigetdir('Select directory.');
    
    path_out = [pname '_cropped'];
    mkdir(path_out)
    
    files = dir([pname filesep '*.TIF']);
    
    % filter stuff
    r_filter = 30;
    f_filter = fspecial('gaussian', 2*r_filter , r_filter); % gaussian filter

    
    for i=1:length(files)
        
        tmp = imread([pname filesep files(i).name]);
        if filter_bool
            tmp2 =  double(tmp(1:2048, 1:2048))-double(imfilter(tmp(1:2048, 1:2048), f_filter, 'same'));
            tmp2 = tmp2-min(tmp2(:));
            imwrite(uint16( tmp2*(2^16-1)./max(tmp2(:)) ), [path_out filesep files(i).name]);
        else
            imwrite(tmp(1:2048, 1:2048), [path_out filesep files(i).name]);
        end
    
    end
   
    disp('Micrographs cropped.')

end
