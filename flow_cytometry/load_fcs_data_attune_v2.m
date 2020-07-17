function [out] = load_fcs_data_attune_v2(pname, fname, radius)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[fcsdat, fcshdr, fcsdatscaled, fcsdat_comp]  = fca_readfcs([pname fname]);

out.fcsdat = fcsdat;
out.fcshdr = fcshdr;
%out.fcsdatscaled = fcsdatscaled;
%out.fcsdat_comp = fcsdat_comp;


%% calculate 
i=2; j=3; %FSC-A vs SSC-A

xy = real([log10(fcsdat(:,i)), log(fcsdat(:,j))]);
NN = get_NN_density_fast(xy, radius);

out.NN_radius = radius;
out.NN = NN;




end

