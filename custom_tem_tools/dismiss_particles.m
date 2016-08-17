function [  ] = dismiss_particles( filter_bool )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    [fname, pname] = uigetfile('*.img', 'Select img stack');

    img = ReadImagic([pname fname]);
    N_particles = size(img,3);

    use = zeros(N_particles, 1);
    r_filter = 8;
    f_filter = fspecial('average', r_filter); % gaussian filter
    
    cf = figure(1);
    go_on = 1;
    i = 1;
    %key = 1;
    while go_on
        j = i; % consider j-th particle
        cur_img_orig = img(:,:, j);
        if filter_bool
            scale = mean(cur_img_orig(:)) *[1 1] + [+1 -1]* std(cur_img_orig(:));
           % tmp =  double(tmp2)-double(imfilter(tmp2, f_filter, 'same'));
            %tmp = tmp-min(tmp(:));
            %tmp =  tmp*(2^16-1)./max(tmp(:));
            tmp =  double(imfilter(cur_img_orig, f_filter, 'same'));
            
            subplot(1,2,1)
            imagesc(tmp), axis image, colormap gray
            title(['Image ' num2str(i) ' of ' num2str(N_particles)  ])

            subplot(1,2,2)
            imagesc(cur_img_orig), axis image, colormap gray

            
        else
            
            tmp = cur_img_orig;
            
            imagesc(tmp), axis image, colormap gray
            title(['Image ' num2str(i) ' of ' num2str(N_particles)  ])
        end
        
        
        
        
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
    WriteImagic(img(:,:,use==0), [pname fname(1:end-4) '_rest'])
    dlmwrite([pname fname(1:end-4) '_cleaned_history.txt'], [[1:N_particles]' use], '\t')
    save([pname fname(1:end-4) '_cleaned_history.mat'])
    disp([num2str(N_particles-sum(use)) ' of ' num2str(N_particles) ' particles dismissed. ' num2str(sum(use)) ' particles remaining.'])
end

