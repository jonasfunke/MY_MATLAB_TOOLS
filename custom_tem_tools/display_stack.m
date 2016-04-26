function [ index_out ] = display_stack( varargin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    if isempty(varargin)
        
        [fname, pname] = uigetfile('*.mat', 'Select a stack (.mat file)');
        data = load([pname fname]);
        
    else
        data = varargin{1};
    end
    dx = 4;
    [n, p, xhist] = uniform_kernel_density(data.angles(:), 3, 0, 120, 0.1);
    
    [angles_sorted, sort_index] = sortrows(data.angles, +1);
    i = [1:size(data.angles,1)]';

    cf = figure();
    
    subplot(3, 2, 1:4)
    imagesc(data.particles(:,:,sort_index(1))), axis image, colormap gray
   % title(['Ref ' num2str(t) ', particle ' num2str(i) ', cc = ' num2str(data(t).stats(sort_index(1),3))])
    cur_img = gca;
    
    subplot(3, 2, 6)
    ylim = [min(data.angles) max(data.angles)];    
    plot(n, xhist)
    set(gca, 'Ylim', ylim)

    grid on

    subplot(3, 2, 5)
    plot(i, angles_sorted, 'b'), hold on
    set(gca, 'YLim', ylim, 'XLim', [1 i(end)]);
    h = imline(gca, round(length(data.angles)/2).*[1 1], ylim);
    setColor(h,[1 0 0]);
    setPositionConstraintFcn(h, @(pos)[ min( i(end), max(1,[pos(2,1);pos(2,1)])) ylim'   ])

    id = addNewPositionCallback(h, @(pos) update_img(  data.particles(:,:,sort_index(  max(1, min(i(end), round(pos(1,1))))   ) ), cur_img )  );
    id2 = addNewPositionCallback(h, @(pos) title(['Angle = ' num2str( angles_sorted(  max(1, min(i(end), round(pos(1,1))))) )]) );
    grid on

    pos_line = wait(h);
    index_selected = max(1, min(i(end), round(pos_line(1,1))));
    
    index_out = sort_index(index_selected);
    pause(0.01)
    close(cf)




end

