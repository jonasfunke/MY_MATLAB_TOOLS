%%
close all, clear all, clc

%% create data set
data = sqrt(rand(100000,2));
resolution  = 0.05;


%tic
%NN1 = get_NN_density(data, resolution);
%toc

%data(500,1) = Inf;
%data(600,1) = nan;
%data(800,2) = nan;

tic
[NN] = get_NN_density_fast(data, resolution);
toc

figure(3), clf
%subplot(1, 2, 1)
%scatter(data(:,1), data(:,2), 5, NN1)

subplot(1, 2, 2)
scatter(data(:,1), data(:,2), 5, NN)

%% create data set
 data = sqrt(rand(1000,2));
 xy = data;
 xy =[
     2     1
     4     5
     3     1]
 
 
% make a grid
resolution = 0.1;


Nx = ceil((max(xy(:,1))-min(xy(:,1)))/resolution);
Ny = ceil((max(xy(:,2))-min(xy(:,2)))/resolution);

% put data into grid

i = max(ceil( (xy(:,1)-min(xy(:,1)))/resolution),1);
j = max(ceil( (xy(:,2)-min(xy(:,2)))/resolution),1);

xygrid = zeros(Nx, Ny);
for k=1:size(xy,1)
    xygrid(i(k), j(k)) = xygrid(i(k), j(k)) +1;
end
    


%%


% get the NN back
tic
NN = xygrid(sub2ind(size(xygrid), i,j));
toc


%%
figure(2), clf
subplot(2, 2, 1)
scatter(data(:,1), data(:,2), 5, NN)
colorbar
subplot(2, 2, 2)
imagesc(xygrid), colorbar



%%
tic
xygrid = zeros(length(xgrid), length(ygrid));
for k=1:size(xygrid,1)
    for l=1:size(xygrid,2)
        xygrid(k,l) = sum(i==k & j==l);
    end
end
toc
subplot(2, 2, 3)
 scatter(data(:,1), data(:,2), '.')

subplot(2, 2, 4)
imagesc(xygrid), colorbar
