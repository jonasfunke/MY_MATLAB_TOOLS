function [i_gated] = one_region_gate_data(data, i, j, roi_position)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    cf = figure();
    cur_roi = drawpolygon('Position', roi_position,'Color','r');
    i_gated = inROI(cur_roi,data(:,i), data(:,j)) ;
    close(cf)
    pause(0.01)
end

