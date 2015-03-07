%% startup
close all, clear all, clc

path0 = cd;
%%
cd(data_directory)

[fname, pname] = uigetfile('.mat', 'Get file');
data1 = load([pname fname]);

sep = strfind(pname, filesep);
path_out = pname(1:sep(end-1));

cd(path_out);
[fname, pname] = uigetfile('.mat', 'Get file');
data2 = load([pname fname]);

cd(path_out);
[fname, pname] = uigetfile('.mat', 'Get file');
data3 = load([pname fname]);
cd(path0)

%%
path_out = [path_out 'combined_' datestr(now, 'yyyy-mm-dd_HH-MM')];
mkdir(path_out)

%% combine to data struct
yield = [[10:24]' data1.yield];
yield = [yield; [25:39]' data2.yield];
yield = [yield; [40:51]' data3.yield];

%% save data
save([path_out filesep 'data.mat'])

%%
cur_fig = figure();
plot(yield(:,1),  yield(:,2), 'r.-')
set(gca, 'YLim', [0 1])
xlabel('Length of spacer [bp]')
ylabel('Yield of reaction')
legend({'cy5-channel'})
print(cur_fig, '-dtiff','-r500' , [path_out filesep 'Yield_combined.tif']); %save figure


%%
cur_fig = figure();
offset = median(yield(end-5:end-1,2));
y =  yield(:,2)-offset;
y = y./ median(y(1:10));
plot(yield(:,1), y, 'r.-')
set(gca, 'YLim', [0 1.2])
xlabel('Length of spacer [bp]')
ylabel('Yield of reaction')
legend({ 'cy5-channel'})

print(cur_fig, '-dtiff','-r500' , [path_out filesep 'Yield_combined_normalized.tif']); %save figure

%%
disp('Finished')

%% combine to data struct
% yield_10uM = [[0:10, 15, 20]' data1.yield];
% yield_10uM = [yield_10uM; [25, 30, 35, 40, 45, 50, 100]' data2.yield(1:7)];
% 
% yield_2mM = [[0:5]' data2.yield(8:13)];
% yield_2mM = [yield_2mM; [6:10, 15, 20, 25, 30, 35, 40, 45, 50, 100]' data3.yield];

% save([path_out filesep 'data.mat'])

% %%
% i1 = 1;
% i2 = 19;
% cur_fig = figure();
% plot(yield_10uM(i1:i2,1), yield_10uM(i1:i2,2), 'g.-', yield_2mM(i1:i2,1), yield_2mM(i1:i2,2), 'r.-')
% set(gca, 'YLim', [0 0.5])
% xlabel('Length of poly-T leash [bases]')
% ylabel('Yield of reaction')
% legend({'10 uM BMH', '2 mM BMH'})
% 
% print(cur_fig, '-dtiff','-r500' , [path_out filesep 'Yield_combined.tif']); %save figure
% 
% 
% %%
% i1 = 1;
% i2 = 20;
% cur_fig = figure();
% plot(ssDNA_mean_end_to_end_distance(yield_10uM(i1:i2,1), 2.3), yield_10uM(i1:i2,2), 'g.', ...
%     ssDNA_mean_end_to_end_distance( yield_2mM(i1:i2,1), 2.3), yield_2mM(i1:i2,2), 'r.-')
% set(gca, 'YLim', [0 0.5])
% xlabel('Mean length of poly-T leash [bases]')
% ylabel('Yield of reaction')
% legend({'10 uM BMH', '2 mM BMH'})
% 
% print(cur_fig, '-dtiff','-r500' , [path_out filesep 'Yield_combined.tif']); %save figure



%%