%% refine 
disp('Refining particles...')

limit = zeros(n_ref/(mirror+1), 1);
if refine
    for t=1:n_ref/(mirror+1)
        [cc_sort, sort_index] = sortrows(data(t).stats(:,3), -1);
        i = [1:size(cc_sort,1)]';

        close all
        subplot(1, 2, 1)
        imagesc(data(t).particles_rot(:,:,sort_index(1))), axis image, colormap gray
       % title(['Ref ' num2str(t) ', particle ' num2str(i) ', cc = ' num2str(data(t).stats(sort_index(1),3))])
        cur_img = gca;

        subplot(1, 2, 2)
        plot(i, cc_sort, 'b'), hold on
        ylim = [0 1];
        set(gca, 'YLim', ylim, 'XLim', [1 i(end)]);
        h = imline(gca,[1 1], ylim);
        setColor(h,[1 0 0]);
        setPositionConstraintFcn(h, @(pos)[ min( i(end), max(1,[pos(2,1);pos(2,1)])) ylim'   ])

        id = addNewPositionCallback(h, @(pos) update_img(  data(t).particles_rot(:,:,sort_index(  max(1, min(i(end), round(pos(1,1))))   ) ), cur_img )  );
        id2 = addNewPositionCallback(h, @(pos) title(['cc = ' num2str( cc_sort(  max(1, min(i(end), round(pos(1,1))))) )]) );

        pos_line = wait(h);
        limit_index = max(1, min(i(end), round(pos_line(1,1))));
        limit(t) = cc_sort(limit_index);

    end
    pause(0.1)
    close all
end
pause(0.1)
close all
