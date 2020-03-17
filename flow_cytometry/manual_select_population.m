function [roi, i_gated] = manual_select_population(xy)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% Select region of interest (ROI)
%i=2; j=4; %FSC-A vs SSC-A


cur_fig = figure(); clf

NN = get_NN_density_fast(xy, 0.02);

scatter(xy(:,1), xy(:,2) ,5, NN, '.'), hold on
colorbar
xlabel('X'), ylabel('Y')
set(gca,'xscale','log','yscale','log')
grid on

title('Select Region of Interest')
roi = drawpolygon;
i_gated = inROI(roi, xy(:,1), xy(:,2) ) ;
pause(0.1)
close(cur_fig)




end

