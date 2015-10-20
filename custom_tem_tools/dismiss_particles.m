function [  ] = dismiss_particles(  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    [fname, pname] = uigetfile('*.img', 'Select img stack');

    img = ReadImagic([pname fname]);
    N_particles = size(img,3);

    use = zeros(N_particles, 1);

    cf = figure(1);
    go_on = 1;
    i = 1;
    %key = 1;
    while go_on
        j = i; % consider j-th particle

        imagesc( img(:,:, j)), axis image, colormap gray
        title(['Image ' num2str(i) ' of ' num2str(N_particles)  ])

        k=waitforbuttonpress;
        tmp = get(gcf,'currentcharacter');
        if strcmp(tmp, 'd')
            use(j) = 0;
            disp([num2str(j) ' dismissed'])    
        else
            if strcmp(tmp, 'b')
                i = i-2;
                disp('back')
            else
                if strcmp(tmp, 'q')
                     go_on = 0;
                     disp('Quit')
                else
                    if strcmp(tmp, 'w')
                        save([pname fname(1:end-4) '_sorting_1-' num2str(i) '.mat'], 'use')
                        disp('data written')
                    else
                        disp([num2str(j) ' used'])
                        use(j) = 1;
                    end
                end
            end
        end
        pause(0.01)

        if i==N_particles
            go_on = 0;
        else
            i=i+1;
        end
    end

    close(cf)
    % write output
    WriteImagic(img(:,:,use==1), [pname fname(1:end-4) '_cleaned'])
    dlmwrite([pname fname(1:end-4) '_cleaned_history.txt'], [[1:N_particles]' use], '\t')
end
