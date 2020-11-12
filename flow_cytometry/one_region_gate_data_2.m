function [i_gated] = one_region_gate_data_2(xy, roi_position)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    cf = figure();
    %scatter(xy(:,1), xy(:,2) ,5, '.'), hold on
    cur_roi = drawpolygon('Position', roi_position,'Color','r');
    i_gated = inROI(cur_roi,xy(:,1), xy(:,2)) ;
    %pause
    close(cf)
    pause(0.01)
end

