function [ data_out, data_out2 ] = sort_stack( varargin )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    if isempty(varargin)
        [fname, pname] = uigetfile('*.mat', 'Select a stack (.mat file)');
        data = load([pname fname]);
    else
        data = varargin{1};
      %  pname = cd;
    end
    
    N_particles = size(data.particles,3);
    use = zeros(N_particles, 1);
    
    cf = figure(1);
    go_on = 1;
    i = 1;
    while go_on
        j = i; % consider j-th particle

        imagesc(data.particles(:,:, j)), axis image, colormap gray
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
                   % if strcmp(tmp, 'w')
                        %save([pname fname(1:end-4) '_sorting_1-' num2str(i) '.mat'], 'use')
                        %disp('data written')
                   % else
                        disp([num2str(j) ' used'])
                        use(j) = 1;
                   % end
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
    
    data_out = data; % copy struct
    data_out.angles = data.angles(use==1);
    data_out.particles = data.particles(:,:,use==1);
    if isempty(data.history)
        data_out.history = {['sorted on ' datestr(now, 'yyyy-mm-dd_HH-MM')]};
    else
        data_out.history = [data_out.history; {['sorted on ' datestr(now, 'yyyy-mm-dd_HH-MM')]}];
    end
    message = [num2str(N_particles-sum(use)) ' of ' num2str(N_particles) ' particles dismissed. ' num2str(sum(use)) ' particles remain.'];
    
    data_out.history = [data_out.history; {message}];
    disp(message)

    
    data_out2 = data; % copy struct
    data_out2.angles = data.angles(use==0);
    data_out2.particles = data.particles(:,:,use==0);
    if isempty(data.history)
        data_out2.history = {['sorted on ' datestr(now, 'yyyy-mm-dd_HH-MM')]};
    else
        data_out2.history = [data_out2.history; {['sorted on ' datestr(now, 'yyyy-mm-dd_HH-MM')]}];
    end
    message = [num2str(N_particles-sum(use)) ' of ' num2str(N_particles) ' particles dismissed. Dismissed particles remain in this stack. ' num2str(sum(use)) ' particles.'];
    
    data_out2.history = [data_out2.history; {message}];
    disp(message)
    
    
    close(cf)
    

end

