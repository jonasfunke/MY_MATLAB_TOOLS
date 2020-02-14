function [] = make_nice_plot(figure_size)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

legend boxoff 
set(gca,'box','off', 'TickDir','out')
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 figure_size(1) figure_size(2) ], 'PaperSize', figure_size );

end

