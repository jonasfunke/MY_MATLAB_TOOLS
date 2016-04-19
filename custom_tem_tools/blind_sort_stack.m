function [ data_out ] = blind_sort_stack( data_array, output_file, filter_bool )
%
%   Detailed explanation goes here


    % compine data sets
    all_angles = zeros(0,3);
    for i=1:length(data_array)
        all_angles = [all_angles; data_array(i).angles, [1:length(data_array(i).angles)]', i.*ones(length(data_array(i).angles),1)];
    end
    
    %
    N_particles = size(all_angles,1);
    p = randperm(N_particles)';
    use = zeros(N_particles, 1);
    
    
    if filter_bool
       r_filter = 20;
       f_filter = fspecial('gaussian', 2*r_filter , r_filter); % gaussian filter
    end
    % make a waitbar
    h = waitbar(0,'Completed');
   
    % loop through mixed data set
    cf = figure(1);
    go_on = 1;
    i = 1;
    while go_on
        j = p(i); % consider j-th particle

        if filter_bool
            img = data_array(all_angles(j,3)).particles(:,:,all_angles(j,2)); 

            %subplot(1, 2, 1)
            %imagesc( double(img)-double(imfilter(img, f_filter, 'same')) ), axis image, colormap gray, colorbar
            
            %subplot(1, 2, 2)
            img_filt = double(img)-double(imfilter(img, f_filter, 'same'));
            imagesc( img_filt , mean(img_filt(:))*[1 1]+[-4 +4]*std(img_filt(:))), axis image, colormap gray

            
        else

             img = data_array(all_angles(j,3)).particles(:,:,all_angles(j,2));
             imagesc(img, mean(img(:))*[1 1]+[-4 +4]*std(img(:))  ), axis image, colormap gray
        end
        title(['Image ' num2str(i) ' of ' num2str(N_particles) ', angle = ' num2str(all_angles(j,1)) ])

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
                        save([output_file '_sorting_1-' num2str(i) '.mat'], 'use', 'p')
                        disp('data written')
                        i = i-1;
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
        waitbar(i/N_particles, h, [num2str(round(100*i/N_particles)) '% completed'])
    end
    
    used_angles = all_angles(use==1,:);    
    data_out = data_array; % copy struct
    for i=1:length(data_array)
        data_out(i).angles = used_angles(used_angles(:,3)==i,1);
        data_out(i).particles =  data_array(i).particles(:,:, used_angles(used_angles(:,3)==i, 2));
   
        if isempty(data_array(i).history)
            data_out(i).history = {['sorted on ' datestr(now, 'yyyy-mm-dd_HH-MM')]};
        else
            data_out(i).history = [data_out(i).history; {['sorted on ' datestr(now, 'yyyy-mm-dd_HH-MM')]}];
        end
        message = [num2str(length(data_array(i).angles)-length(data_out(i).angles)) ...
            ' of ' num2str(length(data_array(i).angles)) ' particles dismissed. ' num2str(length(data_out(i).angles)) ' particles remain.'];

        data_out(i).history = [data_out(i).history; {message}];
        disp(message)
    end
    save([output_file '_sorting_all.mat'], 'data_out', 'use', 'p')
    
    close(h)
    close(cf)
    

end

