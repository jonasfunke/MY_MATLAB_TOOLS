function [  ] = filter_tem_micrographs( r_filter )
%crop micrographs

    [ pname] =uigetdir('Select directory.');
    
    path_out = [pname '_filtered'];
    mkdir(path_out)
    
    files = dir([pname filesep '*.TIF']);
    
    % filter stuff
    f_filter = fspecial('gaussian', 2*r_filter , r_filter); % gaussian filter

    
    for i=1:length(files)
        
        tmp = imread([pname filesep files(i).name]);
        tmp2 =  double(tmp)-double(imfilter(tmp, f_filter, 'same'));
        tmp2 = tmp2-min(tmp2(:));
        imwrite(uint16( tmp2*(2^16-1)./max(tmp2(:)) ), [path_out filesep files(i).name]);

    end
   
    disp([num2str(length(files)) ' micrographs filtered.'])

end
