function [i_selected, roi_pos] = one_region_gate(data, i, j, radius)
%Gate based on columns i and j

cur_fig = figure(); clf

NN = get_NN_density_fast(data(:,[i j]), radius);

scatter(data(:,i), data(:,j), 5, NN, '.'), hold on
colorbar
xlabel('X'), ylabel('Y')
set(gca,'xscale','log','yscale','log')
grid on

title('Select points')
roi1 = drawpolygon;
i_selected = inROI(roi1, data(:,i), data(:,j) ) ;
roi_pos =  roi1.Position;

pause(0.1)
close(cur_fig)




end

