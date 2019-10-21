function [out] = load_fcs_data(pname, fname, path_out, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[fcsdat, fcshdr, fcsdatscaled, fcsdat_comp]  = fca_readfcs([pname fname]);

out.fcsdat = fcsdat;
out.fcshdr = fcshdr;
out.fcsdatscaled = fcsdatscaled;
out.fcsdat_comp = fcsdat_comp;

%% Select region of interest (ROI)
i=2; j=4; %FSC-A vs SSC-A


cur_fig = figure(); clf
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 15 15 ], 'PaperSize', [15 15] );

scatter(fcsdat(:,i), fcsdat(:,j), 1, '.'), hold on
xlabel(fcshdr.par(i).name), ylabel(fcshdr.par(j).name)
set(gca,'xscale','log','yscale','log')
grid on

if  isempty(varargin)
    title('Select Region of Interest')
    roi = drawpolygon;
    i_gated = inROI(roi,fcsdat(:,i), fcsdat(:,j)) ;
    scatter(fcsdat(i_gated,i), fcsdat(i_gated,j),1, 'r.') 
    out.roi_position = roi.Position;

else
    cur_roi = drawpolygon('Position', varargin{1});
    
    i_gated = inROI(cur_roi,fcsdat(:,i), fcsdat(:,j)) ;
    scatter(fcsdat(i_gated,i), fcsdat(i_gated,j),1, 'r.') 
    out.roi_position = varargin{1};
end
out.i_gated = i_gated;

title(fname(1:end-4))

xlabel(fcshdr.par(i).name), ylabel(fcshdr.par(j).name)
legend({'Raw', 'Gated'}, 'location', 'best')

print(cur_fig, '-dpdf', [path_out fname(1:end-4) '_gated.pdf']); %save figure

pause(0.1)

close(cur_fig)




end

