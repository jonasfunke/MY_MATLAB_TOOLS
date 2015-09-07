function [  ] = combine_img_files(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [ pname] = uigetdir('Select directory.');
    
    path_out = pname;

    files = dir([pname filesep '*.img']);
    tmp = ReadImagic([pname filesep files(1).name]);
    
    % compute number of paticles
    N_particles = 0;
    for i=1:length(files)
       info = ReadImagicHeader([pname filesep files(i).name]);
       N_particles = N_particles + size(info.Vals,2);
       tmp= ReadImagic([pname filesep files(i).name]);    
    end
    disp(['You have ' num2str(N_particles) ' particles'])
    
    img = zeros(size(tmp,1),size(tmp,2), N_particles, class(tmp));
	% read images
    counter = 1;
    for i=1:length(files)
         tmp= ReadImagic([pname filesep files(i).name]);    
         for j=1:size(tmp,3)
            img(:,:,counter) = tmp(:,:,j);
            counter = counter+1;
         end
    end
    %% write img file
    WriteImagic(img, [pname filesep 'combined_particles'])
    disp('Particles written to file.')

end

