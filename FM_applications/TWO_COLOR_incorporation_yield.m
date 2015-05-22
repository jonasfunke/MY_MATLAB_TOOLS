%% startup
clc, clear all, close all
path0 = cd;
run('my_prefs.m')


%% choose colors
rgb={'red','green','blue'};
[colors,ok]=listdlg('PromptString', 'Select two colors to be analyzed',...
                'ListString', rgb,...
                'OKString', 'Engage');
while ne(length(colors),2)
    [colors,ok]=listdlg('PromptString', 'Select _TWO_ colors to be analyzed',...
                'ListString', rgb,...
                'OKString', 'Engage');
end

channel = cell(2,1);
channel{1} = rgb{colors(1)};
channel{2} = rgb{colors(2)};

%% LOAD STACK OF MOVIES
pname=uigetdir(data_dir,'Choose the folder with all .fits files.');
files_ch1 = pickFirstFitsFiles(pname, channel{1}); 
files_ch2 = pickFirstFitsFiles(pname, channel{2});

N_movie = length(files_ch1);
if size(files_ch1,1) ~= size(files_ch2,1)
    disp('WARNING: not same number of movie files!')
end

path_out = [pname filesep datestr(now, 'yyyy-mm-dd_HH-MM') '_analysis'];
mkdir(path_out)



%%
%% SET PARAMETER
input = {'First Frame:', 'Last Frame (-1=all):', ['Sequence ' channel{1} ':'], ['Sequence ' channel{2} ':'],... % sample options
    'Radius of peak [pixel]:', 'Integration radius [pixel]:', 'Frames for average [frames]:'};
input_default = {'2', '-1', '01', '10', '4', '3', '-1'};
tmp = inputdlg(input, 'Parameters', 1, input_default);

first = round(str2double(tmp(1))); % first image to read from file
last = round(str2double(tmp(2))); % last image to read from file
%determine sequences 
sequence_ch1 = zeros(1, size(tmp{3},2));
for i=1:size(tmp{3},2)
    if(tmp{3}(i) == '1')
        sequence_ch1(1,i) =1;
    end
end
sequence_ch2 = zeros(1, size(tmp{4},2));
for i=1:size(tmp{4},2)
    if(tmp{4}(i) == '1')
        sequence_ch2(1,i) =1;
    end
end
r_find = str2double(tmp(5)); % radius used to find spots
r_integrate = str2double(tmp(6)); % radius used for integration of intesituies
N_frames = str2double(tmp(7)); % minimal number of found spots in a trace


%% generate movie classes
ch1 = cell(N_movie,1);
ch2 = cell(N_movie,1);

for i=1:N_movie
    ch1{i} = movie(pname, files_ch1{i}, first, last, sequence_ch1); % pname, fname, first, last, sequence
    ch2{i} = movie(pname, files_ch2{i}, first, last, sequence_ch2); % pname, fname, first, last, sequence
end%% generate movie classes

%% make average images

avg_img = cell(N_movie, 2);

for i=1:N_movie
    avg_img{i, 1} = ch1{i}.average_image(N_frames);
    avg_img{i, 2} = ch2{i}.average_image(N_frames);
end


%% determine thresholds and find peaks
peaks_raw = zeros(0,5);

for i=1:N_movie
    [h_min, pch1] = ch1{i}.get_h_min(r_find, N_frames);
    [h_min, pch2] = ch2{i}.get_h_min(r_find, N_frames);
    
    % map peaks
    trace_map = map_traces(pch1(:,1:2), pch2(:,1:2), pch2(:,1:2), r_find*2)+1; %map the tarces from averga positions

    tmp = zeros(size(trace_map,1),5);
    % combine pairs
    for j=1:size(trace_map,1)
        tmp(j,:) = [pch1(trace_map(j,1), 1:2)+1 pch2(trace_map(j,2), 1:2)+1 i]; %x_1 y_1 x_2 y_2 frame
    end
    
    peaks_raw = [peaks_raw; tmp];
end

N_peaks_raw = size(peaks_raw,1);
display(['You have ' num2str(N_peaks_raw) ' pairs'])

%% trace movies

for i=1:N_movie
    bla1 = ch1{1}.traces_movie_position(peaks_raw(:,1:2), 2);
    bla2 = ch2{1}.traces_movie_position(peaks_raw(:,3:4), 2);
   
    for j=1:N_peaks_raw
        plot(bla1{j}(:,1), bla1{j}(:,4), 'r', bla2{j}(:,1), bla2{j}(:,4), 'g' )
        pause
    end
end




%% determinet find peaks in subframes
i=1;




img1 = avg_img{i,1}(128:128+256-1, 128:128+256-1);      
p1 = find_peaks2d(img1, r_find, 0, 0); % finding all possible peaks p has x, y, height, height-bg, I, I-I_bg
p1(:,1:2) = p1(:,1:2)+1;
img1_mean = mean(img1(:));
img1_std = std(img1(:)); 

img2 = avg_img{i,2}(128:128+256-1, 128:128+256-1);          
p2 = find_peaks2d(img2, r_find, 0, 0); % finding all possible peaks p has x, y, height, height-bg, I, I-I_bg
p2(:,1:2) = p2(:,1:2)+1;
img2_mean = mean(img2(:));
img2_std = std(img2(:)); 


close all
figure('units','normalized','outerposition',[0 0 1 1])





subplot(3,4, [1 2 5 6])
plot_image(img1, [0.1 0.3]), hold on
cur_img1 = gca;


subplot(3,4,[9 10])
xhist = min(p1(:,4)):5:max(p1(:,4));
n = hist(p1(:,4), xhist);
semilogy(xhist, sum(n)-cumsum(n)), hold on

ylim = [1 size(p1,1)];
xlim = [min(p1(:,4)) max(p1(:,4))];
set(gca, 'YLim', ylim, 'XLim', xlim)
plot(p1(p1(:,4)>=img_mean, 1), p1(p1(:,4)>=img_mean, 2), 'ro', 'parent', cur_img1)
h = imline(gca,[img_mean img_mean], ylim);
setColor(h,[1 0 0]);
setPositionConstraintFcn(h, @(pos) [ min( xlim(2), max(1,[pos(2,1);pos(2,1)])) ylim'   ])
id = addNewPositionCallback(h, @(pos) update_peaks( p1, pos(1,1)  , cur_img1, 'ro' )  ); % plot circles 
%pos_line = wait(h);
%h_min = pos_line(1,1); % minimal height of peak

subplot(3,4,[3 4 7 8])
plot_image(img2, [0.1 0.3]), hold on
cur_img2 = gca;

subplot(3,4,[11 12])
xhist = min(p2(:,4)):5:max(p2(:,4));
n = hist(p2(:,4), xhist);
semilogy(xhist, sum(n)-cumsum(n)), hold on

ylim = [1 size(p2,1)];
xlim = [min(p2(:,4)) max(p2(:,4))];
set(gca, 'YLim', ylim, 'XLim', xlim)
plot(p2(p2(:,4)>=img_mean, 1), p2(p2(:,4)>=img_mean, 2), 'ro', 'parent', cur_img2)
h = imline(gca,[img_mean img_mean], ylim);
setColor(h,[0 1 0]);
setPositionConstraintFcn(h, @(pos) [ min( xlim(2), max(1,[pos(2,1);pos(2,1)])) ylim'   ])
id = addNewPositionCallback(h, @(pos) update_peaks( p2, pos(1,1)  , cur_img2, 'go') ); % plot circles 
id2 = addNewPositionCallback(h, @(pos) update_peaks( p2, pos(1,1)  , cur_img1, 'go') ); % plot circles 

pos_line = wait(h);
h_min = pos_line(1,1); % minimal height of peak


p_out = p(p(:,4)>=h_min,:);









