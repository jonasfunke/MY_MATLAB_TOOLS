function [  ] = convert_mrc_stack_to_mrc()
    %filters mrc stack and writes normalized 8bit jpg 
    cur_dir = cd;
    [fname, pname] =uigetfile({'*.mrc'; '*.mrcs'}, 'Select a mrc stack.', cur_dir);
    
    path_out = [pname fname(1:end-4) '_Micrographs'];
    
    mkdir(path_out)
    
    
    
    [images, info] = ReadMRC([pname fname]);
    disp('read')
    for i=1:size(images,3)
        file_out = [path_out filesep fname(1:end-4) sprintf('_%03i.', i)  '.mrc'];
        WriteMRC(single(images(:,:,i)), 1, file_out)
    end
   
    disp(['Micrographs filtered and written to ' path_out])

end
