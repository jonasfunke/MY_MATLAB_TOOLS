function [roi1_pos, roi2_pos, i_gated1, i_gated2] = manual_select_population(xy)
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

title('Select Region 1')
roi1 = drawpolygon;
i_gated1 = inROI(roi1, double(xy(:,1)), double(xy(:,2)) ) ;
roi1_pos =  roi1.Position;

% dead_pgon = polyshape(roi1_pos);
% plot(dead_pgon, 'FaceColor', 'none', 'EdgeColor', 'r')

title('Select Region 2')
roi2 = drawpolygon;
i_gated2 = inROI(roi2, double(xy(:,1)), double(xy(:,2)) ) ;
roi2_pos =  roi2.Position;


pause(0.1)
close(cur_fig)




end

