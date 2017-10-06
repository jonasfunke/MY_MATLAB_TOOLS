function [  ] = filter_mrc_stack()
    %filters mrc stack and writes normalized 8bit jpg 
    cur_dir = cd;
    [fname, pname] =uigetfile({'*.mrc'; '*.mrcs'}, 'Select a mrc stack.', cur_dir);
    
    path_out = [pname fname(1:end-4) '_filtered'];
    
    mkdir(path_out)
    
    
    % filter stuff
    r_filter = 50;
    f_filter = fspecial('gaussian', 2*r_filter , r_filter); % gaussian filter
   
    format_out = 'jpg';
    bit_depth = 8;
    
    disp(['Filtering with high-pass gaussian radius of ' num2str(r_filter) ' px'])
    
    disp(['Output: ' num2str(bit_depth) ' bit ' num2str(format_out)])
    
    [images, info] = ReadMRC([pname fname]);
    
    for i=1:size(images,3)
        file_out = [path_out filesep fname(1:end-4) sprintf('_%03i.', i)  format_out];
        tmp =  double(images(:,:,i))-double(imfilter(images(:,:,i), f_filter, 'same', 'replicate'));
        tmp = tmp-min(tmp(:));
        imwrite(uint8( tmp*(2^bit_depth-1)./max(tmp(:)) ), file_out, format_out,'BitDepth',bit_depth);
    end
   
    disp(['Micrographs filtered and written to ' path_out])

end
