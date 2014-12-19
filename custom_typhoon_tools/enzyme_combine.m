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
path_out = [path_out 'combined'];
mkdir(path_out)
%% combine to data struct
yield = [[10:24]' data1.yield];
yield = [yield; [25:39]' data2.yield];
yield = [yield; [40:51]' data3.yield];
%%
cur_fig = figure();
plot(yield(:,1), yield(:,2), 'g.-', yield(:,1), yield(:,3), 'r.-')
set(gca, 'YLim', [0 0.5])
xlabel('Length of spacer [bp]')
ylabel('Yield of reaction')
legend({'cy3-channel', 'cy5-channel'})

print(cur_fig, '-dtiff','-r500' , [path_out filesep 'Yield_combined.tif']); %save figure


%%
cur_fig = figure();
offset = median(yield(end-5:end-1,3));
y =  yield(:,3)-offset;
y = y./ median(y(1:10));
plot(yield(:,1), y, 'r.-')
set(gca, 'YLim', [0 1.2])
xlabel('Length of spacer [bp]')
ylabel('Yield of reaction')
legend({'cy3-channel', 'cy5-channel'})

print(cur_fig, '-dtiff','-r500' , [path_out filesep 'Yield_combined.tif']); %save figure


%%