function [  ] = combine_mrcs_files( pixA )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
    [ pname] = uigetdir('Select directory.');
    
    path_out = pname;

    files = dir([pname filesep '*.mrcs']);
    tmp = ReadMRC([pname filesep files(1).name]);
  
   
%%
    img = zeros(size(tmp,1),size(tmp,2), 0, 'single');
	
    % read images
    for i=1:length(files)
         tmp = ReadMRC([pname filesep files(i).name]);    
         img = cat(3, img, tmp);
         disp(['Reading stack ' num2str(i) ' of ' num2str(length(files))])
    end
    disp([num2str(size(img,3)) ' particles read.'])
    
    %% write img file
    WriteMRC(img, pixA, [pname filesep 'combined_particles.mrcs'], 2, size(img,3) )
    disp(['Particles written to ' pname filesep 'combined_particles.mrc'])
%%
end