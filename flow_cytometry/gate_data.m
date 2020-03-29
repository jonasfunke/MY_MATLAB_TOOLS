function [i_gated] = gate_data(xy, roi_position)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    cf = figure();
    cur_roi = drawpolygon('Position', roi_position,'Color','r');
    i_gated = inROI(cur_roi,xy(:,1), xy(:,2)) ;
    close(cf)
    pause(0.01)
end

