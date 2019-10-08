function [out] = load_fcs_data(pname, fname, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[fcsdat, fcshdr, fcsdatscaled, fcsdat_comp]  = fca_readfcs([pname fname]);

out.fcsdat = fcsdat;
out.fcshdr = fcshdr;
out.fcsdatscaled = fcsdatscaled;
out.fcsdat_comp = fcsdat_comp;

%% Select region of interest (ROI)
cur_fig = figure(2); clf

i=2; j=4; %FSC-A vs SSC-A
scatter(fcsdat(:,i), fcsdat(:,j), 1, '.'), hold on
xlabel(fcshdr.par(i).name), ylabel(fcshdr.par(j).name)
set(gca,'xscale','log','yscale','log')
grid on
title('Select Region of Interest')

if  length(varargin)==0
    roi = drawpolygon;
    i_gated = inROI(roi,fcsdat(:,i), fcsdat(:,j)) ;
    scatter(fcsdat(i_gated,i), fcsdat(i_gated,j), 'r.') 
    out.roi = roi;
else
    i_gated = inROI(roi,fcsdat(:,i), fcsdat(:,j)) ;
    scatter(fcsdat(i_gated,i), fcsdat(i_gated,j), 'r.') 
    out.roi = varargin{1};
end

pause
close(2)


out.i_gated = i_gated;

end

