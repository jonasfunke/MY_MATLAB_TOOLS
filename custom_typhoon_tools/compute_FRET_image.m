function [  ] = compute_FRET_image(floc, scale, normalization, size, path_out )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

   % [fname, pname] = uigetfile('*.mat', 'Select mat file');
    %floc = [pname fname];
    load(floc, 'gelData', 'gamma_calc');
    
    DD = gelData.images{1};
    DA = gelData.images{3};
    AA = gelData.images{2};
   % img_sz = size(DD);

    E = DA./(DA+gamma_calc.*DD);

    
        %area = [1610, 1875, 1725, 3075];
   % area = [1600, 1900, 1700, 3100];

    cf = figure;
    set(gcf,'Visible','off', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
        'PaperPosition', [0 0 size(1) size(2)]);


    fh = imagesc(E,  scale); colorbar, axis image
    colormap(flipud(parula))
    opa = min(AA./normalization,1);
    set(fh, 'AlphaData', opa)
    %set(gca, 'Ylim', [area(1) area(2)],  'Xlim', [area(3) area(4)], ...
    %    'Xtick', [], 'YTick', [])

    print(cf, [path_out  gelData.filenames{1}(1:end-4) '_FRET_image_' num2str(round(normalization)) '_' num2str(scale(1)) '-' num2str(scale(2)) '.png'], '-dpng','-r1000' )

    
    
     cf = figure;
    set(gcf,'Visible','off', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
        'PaperPosition', [0 0 size(1) size(2)]);

    imagesc(opa,  [0 1]); colorbar, axis image
    colormap(flipud(parula))

    print(cf, [path_out  gelData.filenames{1}(1:end-4) '_FRET_image_' num2str(round(normalization)) '_' num2str(scale(1)) '-' num2str(scale(2)) '_opacityMask.png'], '-dpng','-r300' )

    disp('done')

end

