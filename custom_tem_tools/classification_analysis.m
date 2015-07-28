%%
close all, clear all, clc

load('/Users/jonasfunke/Dropbox/POSITIONING/Figure_1/angle_calibration_02/data/average_angles.mat')
load('/Volumes/matlabuser/jonasfunke/data/NoInt_FoB25/2015-07-27_14-00_automatic_classification/data_images_all.mat')

no_int = load('/Users/jonasfunke/Dropbox/Nucleosome_Dropbox/Figure4/TEM/Angle_references/no_interaction/iter940/Mean_angles_no_int.txt');
%%
tmp = angle(1:2:end,1);
xhist = 1:38;
n = hist(classification(:,2), xhist);

figure(1)
bar(xhist, n)

figure(2)
alpha = [angle(1:2:end,1); no_int(:,2)];

bar(alpha, n)



