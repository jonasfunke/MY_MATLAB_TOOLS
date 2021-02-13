function [i_gated] = gate_data_2d(xy, roi_position)
% Gate data based on polygon
%   Detailed explanation goes here
    cf = figure();
    %scatter(xy(:,1), xy(:,2) ,5, '.'), hold on
    cur_roi = drawpolygon('Position', roi_position,'Color','r');
    i_gated = inROI(cur_roi, double(xy(:,1)), double(xy(:,2))) ;
    %pause
    close(cf)
    pause(0.01)
end

