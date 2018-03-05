function [  ] = convert_mrc_stack_to_tif()
    %filters mrc stack and writes normalized 8bit jpg 
    cur_dir = cd;
    [fname, pname] =uigetfile({'*.mrcs'; '*.mrc'}, 'Select a mrc stack.', cur_dir);
    
    path_out = [pname fname(1:end-5) '_tif'];
    
    mkdir(path_out)
    
        
    [images, info] = ReadMRC([pname fname]);
    bit_depth = 16;
    format_out = 'tif';
    
    for i=1:size(images,3)
        file_out = [path_out filesep fname(1:end-4) sprintf('_%03i.', i)  format_out];
        tmp = double(images(:,:,i));
        tmp = tmp-min(tmp(:));
        imwrite(uint16( tmp*(2^bit_depth-1)./max(tmp(:)) ), file_out);

        disp([num2str(i) ' of ' num2str(size(images,3)) ' filtered.' ]);
    end
   
    disp(['Micrographs filtered and written to ' path_out])

end
