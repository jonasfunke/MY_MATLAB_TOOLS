function [ leak_dir] = calculate_corrections(dd_bg, da_bg, aa_bg, file_out )
%determines leakage and direct excitation correction from dd_img, da_img
%and aa_img
%   Detailed explanation goes here

%% Integrate D-only and A-only band
%[I, areas] = integrate_areas({da_bg, dd_bg, aa_bg}, 2, 1, [1 0 0]); %cell of images, number of bands, 1=all bands habe the same size
[I, areas] = integrate_areas({dd_bg, da_bg, aa_bg}, 2); %cell of images, number of bands, 


%% make profiles
donly_profiles = [[areas(1,2):areas(1,2)+areas(1,4)]' ...
    sum(dd_bg(areas(1,2):areas(1,2)+areas(1,4), areas(1,1):areas(1,1)+areas(1,3))')' ...
    sum(da_bg(areas(1,2):areas(1,2)+areas(1,4), areas(1,1):areas(1,1)+areas(1,3))')' ...
    sum(aa_bg(areas(1,2):areas(1,2)+areas(1,4), areas(1,1):areas(1,1)+areas(1,3))')'];

aonly_profiles = [[areas(2,2):areas(2,2)+areas(2,4)]' ...
    sum(dd_bg(areas(2,2):areas(2,2)+areas(2,4), areas(2,1):areas(2,1)+areas(2,3))')' ...
    sum(da_bg(areas(2,2):areas(2,2)+areas(2,4), areas(2,1):areas(2,1)+areas(2,3))')' ...
    sum(aa_bg(areas(2,2):areas(2,2)+areas(2,4), areas(2,1):areas(2,1)+areas(2,3))')'];


%% calculate corrections 

%d-only
pos = areas(1,:);
lane_dd = dd_bg(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3));
lane_da = da_bg(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3));

[p_leak, DA_donly, DD_donly] = calculateRation(lane_da, lane_dd, 0);

% a-only
pos = areas(2,:);
lane_aa = aa_bg(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3));
lane_da = da_bg(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3));

[p_dir, DA_aonly, AA_aonly] = calculateRation(lane_da, lane_aa, 0);


%% plot
close all
scrsz = get(0,'ScreenSize');
fig_dim =[2*13.7 2*8];%5.23];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);


subplot(2, 2, 1)
xlim = [min(DD_donly(:)) max(DD_donly(:))];
plot(DD_donly(:), DA_donly(:), 'g.', xlim, p_leak(1)*xlim+p_leak(2), 'k' ,'MarkerSize', 1)
legend({'Data', ['Fit, leak=' num2str(round(p_leak(1)*100)) '%']}, 'FontSize', 10, 'Location', 'SouthEast')
xlabel('D->D')
ylabel('D->A')
title('Leakage correction')
set(gca, 'XLim', [xlim(1) xlim(2)*1.2])

subplot(2, 2, 2)
plot(donly_profiles(:,1), p_leak(1,1)*donly_profiles(:,2), 'g', donly_profiles(:,1), donly_profiles(:,3), 'b', donly_profiles(:,1), donly_profiles(:,3)-donly_profiles(:,2)*p_leak(1,1), 'b--' )
legend({'D->D, scaled', 'D->A', 'D->A corrected'}, 'FontSize', 10)
xlabel('Migration distance [pixel]')
ylabel('Intensity')
set(gca, 'XLim', [donly_profiles(1,1) donly_profiles(end,1)])

subplot(2, 2, 3)
xlim = [min(AA_aonly(:)) max(AA_aonly(:))];
plot(AA_aonly(:), DA_aonly(:), 'r.', xlim, p_dir(1)*xlim+p_dir(2), 'k' ,'MarkerSize', 1)
legend({'Data', ['Fit, dir=' num2str(round(p_dir(1)*100)) '%']}, 'FontSize', 10, 'Location', 'SouthEast')
xlabel('A->A')
ylabel('D->A')
title('Direct excitation correction')
set(gca, 'XLim', [xlim(1) xlim(2)*1.2])

subplot(2, 2, 4)
plot(aonly_profiles(:,1), p_dir(1,1)*aonly_profiles(:,4), 'r', aonly_profiles(:,1), aonly_profiles(:,3), 'b', aonly_profiles(:,1), aonly_profiles(:,3)-aonly_profiles(:,4)*p_dir(1,1), 'b--' )
legend({'A->A, scaled', 'D->A', 'D->A corrected'}, 'FontSize', 10)
xlabel('Migration distance [pixel]')
ylabel('Intensity')
set(gca, 'XLim', [aonly_profiles(1,1) aonly_profiles(end,1)])

print(cur_fig, '-depsc2','-loose' , [file_out(1:end-4) '_corrections']); %save figure

%display('Press some key to continue...')
questdlg('Go on','Halt','Go on','Go on');
close all

%% compine to matrix
leak_dir = [p_leak ; p_dir]; 
disp(['Leakage: ' num2str(round(p_leak(1)*100)) '% Direct Ex.: ' num2str(round(p_dir(1)*100)) '%' ])

%% save data
save(file_out, 'leak_dir' ,'-ascii')

close all


end

