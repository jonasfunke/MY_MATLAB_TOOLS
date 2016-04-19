function [ classes ] = classify_stack(data, class_keys, file_out )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
filter_bool = 1;
    N_particles = length(data.angles);

    classes = zeros(N_particles, 1);
    r_filter = 15;
    f_filter = fspecial('gaussian', 3*r_filter , r_filter); % gaussian filter
    
    cf = figure;
    go_on = 1;
    i = 1;
    
    while go_on
                
        if filter_bool
            tmp2 = data.particles(:,:, i);
            tmp =  double(tmp2)-double(imfilter(tmp2, f_filter, 'same'));
            tmp = tmp-min(tmp(:));
            tmp =  tmp*(2^16-1)./max(tmp(:));
        else
            tmp = data.particles(:,:, i);
        end
        imagesc(tmp), axis image, colormap gray
        title(['Image ' num2str(i) ' of ' num2str(N_particles)  ])

        k=waitforbuttonpress;
        tmp = get(gcf,'currentcharacter');
        
        for j=1:length(class_keys)
            if strcmp(tmp, class_keys{j})
                classes(i) = j;
                disp(['Particle ' num2str(i) ' in class ' class_keys{j}])    
            end
        end
        
        
        if strcmp(tmp, 'q')
             go_on = 0;
             disp('Quit')   
        end
        
        if strcmp(tmp, 'b')
            i = i-2;
            disp('back')
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
    save(file_out, 'data', 'class_keys', 'classes')

    % display
    for i=1:length(class_keys)
        disp(['Class ' class_keys{i} ': ' num2str(sum(classes==i)) ' particles'])
    end
    
end

