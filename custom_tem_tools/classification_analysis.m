%%
close all, clear all, clc


load('/Volumes/matlabuser/jonasfunke/data/NoInt_FoB25/2015-07-27_14-00_automatic_classification/data_images_all.mat')
load('/Users/jonasfunke/Dropbox/POSITIONING/Figure_1/angle_calibration_02/data/average_angles.mat')
no_int = load('/Users/jonasfunke/Dropbox/Nucleosome_Dropbox/Figure4/TEM/Angle_references/no_interaction/iter940/Mean_angles_no_int.txt');



%load('/Volumes/matlabuser/jonasfunke/data/NoInt_FoB25/2015-07-27_14-00_automatic_classification/data_images_all.mat')

% load('/Volumes/matlabuser/jonasfunke/data/NoInt_FoB25/2015-07-28_15-08_automatic_classification/data_before_150531_PK_GR_13__B2_FRETst_noint_stack.mat')
% load('/Volumes/matlabuser/jonasfunke/data/NoInt_FoB25/2015-07-28_15-08_automatic_classification/data_images_all.mat')

load('/Volumes/matlabuser/jonasfunke/data/inner_0deg_Fob25/2015-07-29_09-19_automatic_classification/data_image1-1280.mat')
load('/Volumes/matlabuser/jonasfunke/data/inner_0deg_Fob25/2015-07-29_09-19_automatic_classification/data_before_images_export.mat')


%% correct classification

alpha = 0:dalpha:359;  % angles
n_rot = length(alpha); % number of rotations
classification_cor = zeros(size(classification,1), 5);
cor_matrix_cor = cell(size(cor_matrix));

for i=1:1280 % loop through images

    % determine best rotation
    tmp =min(cor_matrix{i});
    bg = tmp(ones(2*n_rot,1),:);
    cor_matrix_cor{i} = cor_matrix{i}-bg;
    [max_cor, i_tmp] = max(cor_matrix{i}(:)-bg(:)); % use best overall correlation
    [best_rot, best_ref] = ind2sub(size(cor_matrix{i}),i_tmp);
        
    mirror = best_rot > n_rot; % determine if image was mirrored
    % determine rotation angle
    if mirror 
        alpha_found = alpha(best_rot-n_rot);
    else
        alpha_found = alpha(best_rot);
    end
    classification_cor(i,:) = [i best_ref alpha_found mirror max_cor]; % combine output
    
end 


%%
alpha_ref = [angle(1:end,1); no_int(:,2)];

n_hist = zeros(1,N_ref);
for i=1:N_ref
    n_hist(i) = sum(classification_cor(:,2)==i); % count how many images fall into class i
end

%%
figure(1)
plot(n_hist, '.')
figure(2)
bar(alpha_ref, n_hist)
figure(3)
plot(alpha_ref, 1:N_ref, '.')


%%
tmp = [alpha_ref n_hist'];

tmp = sortrows(tmp, 1);

x = tmp(1:end-1,1)+diff(tmp(:,1))/2;

width = diff(x);
%
bla = tmp(2:end-1,2) ./ width;
close all
area(tmp(:,1), tmp(:,2)), hold on
%plot(tmp(2:end-1,1) , bla, '.-')




