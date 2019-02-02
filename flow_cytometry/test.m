%% startup
close all, clear all, clc


%%

floc = '/Users/jonasfunke/Dropbox (Personal)/PlectonicBiotech/Experiments/data_FACS/2018-12-05_JonasFunke_JURKAT_4x105_2_1nM_DO_sdAB.fcs';

%%

[data, param, header] =fcsread(floc);
%%

[fcsdat, fcshdr, fcsdatscaled, fcsdat_comp] = fca_readfcs(floc);
%%
close all
figure(1)
for i=1:7
    subplot(2, 4, i)
    histogram(fcsdat(:,i))
    set(gca,'xscale','log')
    xlabel(fcshdr.par(i).name)
end

%%

%%

loglog(fcsdat(:,1), fcsdat(:,2), '.')

