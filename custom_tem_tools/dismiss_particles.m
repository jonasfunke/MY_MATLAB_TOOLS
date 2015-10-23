function [  ] = dismiss_particles( filter_bool )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    [fname, pname] = uigetfile('*.img', 'Select img stack');

    img = ReadImagic([pname fname]);
    N_particles = size(img,3);

    use = zeros(N_particles, 1);
    r_filter = 15;
    f_filter = fspecial('gaussian', 3*r_filter , r_filter); % gaussian filter
    
    cf = figure(1);
    go_on = 1;
    i = 1;
    %key = 1;
    while go_on
        j = i; % consider j-th particle
        
        if filter_bool
            tmp2 = img(:,:, j);
            tmp =  double(tmp2)-double(imfilter(tmp2, f_filter, 'same'));
            tmp = tmp-min(tmp(:));
            tmp =  tmp*(2^16-1)./max(tmp(:));
        else
            tmp = img(:,:, j);
        end
        imagesc(tmp), axis image, colormap gray
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
    disp([num2str(N_particles-sum(use)) ' of ' num2str(N_particles) ' particles dismissed. ' num2str(sum(use)) ' particles remaining.'])
end

