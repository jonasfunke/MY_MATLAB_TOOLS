%%
close all, clear all, clc

file_loc = '/Users/jonasfunke/Dropbox/Nucleosome_unwrapping/Data/particle_stacks/all_particles/2014-03-12_G5B06_FS-1000+Nucleosome_FoB10+0.0MNaCl_0001-1000_raw_images.mat';
load(file_loc);

colormap gray
%% FFT

img = particles(:,:,1);

r_highpass = 10; % pixel
r_lowpass = 3; % pixel
f_highpass = fspecial('gaussian', 2*r_highpass, r_filter); % gaussian filter
f_lowpass = fspecial('gaussian', 2*r_lowpass, r_filter); % gaussian filter


img_filter =  double(imfilter(img, f_lowpass, 'same'))-double(imfilter(img, f_highpass, 'same'));


img_fourier = abs(fftshift(fft2(img_filter)));

figure(1), clf
subplot(2, 2, 1)
imagesc(img), axis image

subplot(2, 2, 2)
imagesc(img_filter), axis image


subplot(2, 2, 3)
imagesc(abs(fftshift(fft2(img))), [0, 1000]), axis image

subplot(2, 2, 4)
imagesc( img_fourier, [0, 1000]), axis image

%% Cartesian => Polar coordinates


[X,Y] = meshgrid([(-99:100)-0.5], [(-99:100)-0.5]);
%tmp = [X(:), Y(:), img_fourier(:)];

x = X(:); 
y = Y(:);

tmp = img_fourier(:);

t = atan2(X(:),Y(:));
r = sqrt(X(:).^2+Y(:).^2);                  

[T,R] = meshgrid(-pi:0.01:pi,0:200);

A = zeros(size(T));
%%
for i=1:size(A,1)
    for j=1:size(A,2)
        t_cur = T(i,j);
        r_cur = R(i,j);
        
        bla = tmp( (t_cur-0.01 < t) & (t< t_cur+0.01) & (r_cur-1 < r) & (r < r_cur+1));
        if ~isempty(bla)
            A(i,j) = sum(bla)/length(bla);
        else
            A(i,j) = 0;
        end
    end
end

%%
figure(2), clf
subplot(1, 2, 1)
imagesc( img_fourier), axis image

subplot(1, 2, 2)
imagesc(A), axis image
%%
figure(3)
tmp = sum(A);

plot(tmp)
%%
figure(4)
%plot(1:105, tmp(1:105), 1:105, tmp(106:210), 1:105, (tmp(1:105)+tmp(106:210))/2)
plot(1:315, tmp(1:315), 1:315, tmp(315:629) , 1:315, (tmp(1:315)+tmp(315:629))/2)
%%
A=double(img_fourier);   %read image
[m, n]=size(A);



[NN, MM]=meshgrid((-99:100)-0.5,(-99:100)-0.5);
T=atan2(NN,MM);
R=sqrt(MM.^2+NN.^2);                  



[Tq, Rq] = meshgrid(linspace(-pi, pi, 200), linspace(0, 200, 200)); 



Aq = interp2(T, R, A, Tq, Rq,'linear',0);


figure(2), clf
subplot(1, 2, 1)
imagesc( A), axis image

subplot(1, 2, 2)
imagesc(Aq), axis image





%%

A=double(img_fourier);   %read image
[m n]=size(A);
[t r]=meshgrid(linspace(-pi,pi,n),1:m); 

M=m;
N=n;
[NN MM]=meshgrid((1:N)-n/2-0.5,(1:M)-m/2-0.5);
T=atan2(NN,MM);
R=sqrt(MM.^2+NN.^2);                  

B = interp2(t,r,A,T,R,'linear',0);


figure(2), clf
subplot(1, 2, 1)
imagesc( A), axis image

subplot(1, 2, 2)
surf(T,R, B), axis image

%% Average





%% Find peaks




%% 2nd try

%%
close all, clear all, clc

[fname, pname] = uigetfile('*.img', 'Select img stack');

%%
particles = ReadImagic([pname fname]);
N_particles = size(img,3);



%% FFT
for i=10:15
img = particles(:,:,i);

r_highpass =15; % pixel
r_lowpass = 2; % pixel
f_highpass = fspecial('gaussian', 2*r_highpass, r_highpass); % gaussian filter
f_lowpass = fspecial('gaussian', 2*r_lowpass, r_lowpass); % gaussian filter


img_filter =  double(imfilter(img, f_lowpass, 'same'))-double(imfilter(img, f_highpass, 'same'));


img_fourier = abs(fftshift(fft2(img_filter)));

figure(1), clf
subplot(2, 2, 1)
imagesc(img), axis image

subplot(2, 2, 2)
imagesc(img_filter), axis image

subplot(2, 2, 3)
imagesc(abs(fftshift(fft2(img))), [0, 1000]), axis image

subplot(2, 2, 4)
imagesc( img_fourier, [0, 1000]), axis image


[X,Y] = meshgrid([-99.5:99.5]', [-99.5:99.5]');


xy = [X(:), Y(:), img_fourier(:)];

rt = [sqrt(xy(:,1).^2+xy(:,2).^2), atan2(xy(:,1), xy(:,2)), xy(:,3)];

r_max = 50;
rt_2 = rt(rt(:,1)<40,:);


%rt_2 = rt(rt(:,3)>400,:);


%rt_2(rt_2(:,2)<0,2) = rt_2(rt_2(:,2)<0,2)+pi; 

x_points = -pi:0.0175:pi;
n = zeros(size(x_points));
h = 0.0175*2;
for j=1:length(rt_2(:,2)) % loop over data
    n = n + rt_2(j,3)*normpdf(x_points, rt_2(j,2), h);
end

p = n/200; %/sum(rt_2(:,3));



figure(2)
plot(rt_2(:,2)*180/pi,rt_2(:,3), '.', x_points*180/pi, p)



pause
end



