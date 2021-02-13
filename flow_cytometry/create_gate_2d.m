function [roi1_pos] = create_gate_2d(xy, r, axis_labels, command)
% create a gate for points in xy
    %Select region of interest (ROI)

    %xy = log10(xy(xy(:,1)>0 & xy(:,2)>0,:));
    cur_fig = figure(); clf
    
    NN = get_NN_density_fast(real(log10(xy)), r); % calculate neearest neighour data

    scatter(xy(:,1), xy(:,2), 5, NN, '.'), hold on
    colorbar
    xlabel(axis_labels{1}), ylabel(axis_labels{2})
    set(gca,'xscale','log','yscale','log')
    grid on

    title(command)
    roi1 = drawpolygon;
    %i_gated1 = inROI(roi1, xy(:,1), xy(:,2) ) ;
    roi1_pos =  roi1.Position;



    pause(0.1)
    close(cur_fig)




end

