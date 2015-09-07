function [  ] = dismiss_particles(  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    [fname, pname] = uigetfile('*.img', 'Select img stack');
    path_out = pname;

    img = ReadImagic([pname fname]);
    N_particles = size(img,3);

    use = zeros(N_particles, 1);

    figure(1)
    go_on = 1;
    j = 1;
    while go_on

        imagesc( img(:,:,j) ), axis image, colormap gray
        title(['Image ' num2str(j) ' of ' num2str(N_particles) ])

        k = getkey();

        if isempty(k)
            go_on = 0;
        else
            if k==100 % this is d=delete/dismiss
                use(j) = 0;
                disp([num2str(j) ' dismissed'])
            else
                if k==28 % this is back button (left arrow)
                    j = j-2;
                    disp('back')
                else
                    disp([num2str(j) ' accepted'])
                    use(j) = 1;
                end
            end
        end

        if j==N_particles
            go_on = 0;
        else
            j=j+1;
        end
    end
    
    % write output
    WriteImagic(img(:,:,use==1), [pname fname(1:end-4) '_cleaned'])
    dlmwrite([pname fname(1:end-4) '_cleaned_history.txt'], [[1:N_particles]' use], '\t')
end

