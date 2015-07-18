%% description
%% startup
close all, clear all, clc
dalpha = 10;

%% load stack of particles
[fname, pname] = uigetfile('*.img', 'Select the stack of images to be classified');
img = ReadImagic([pname, fname(1:end-4)]); % read imagic file

%% load reference classes
[fname_ref, pname_ref] = uigetfile('*.img', 'Select the stack of reference classes.');
ref = ReadImagic([pname_ref, fname_ref(1:end-4)]); % read imagic file

N_img = size(img,3);
N_ref = size(ref,3);
N = size(img,1);
%% filter references

r_filter = 20;
f_filter = fspecial('gaussian', r_filter , r_filter); % gaussian filter
ref_filtered = zeros(size(ref,1),size(ref,2),N_ref,'single');
subt  = zeros(size(ref,1),size(ref,2),1,N_ref); %, 'single');
for r=1:N_ref
    ref_filtered(:,:,r) = double(ref(:,:,r))-double(imfilter(ref(:,:,r), f_filter, 'same'));
    subt(:,:,1,r) = double(imfilter(ref(:,:,r), f_filter, 'same'));
end

% plot references
bla = zeros(size(ref,1),size(ref,2),1,N_ref, 'single');
for r=1:N_ref
    bla(:,:,1,r) = ref(:,:,r);
end
%figure(1)
%montage(bla, 'DisplayRange', [min(bla(:)) max(bla(:))])

figure(2)
%montage(ref_filtered, 'DisplayRange', [min(ref_filtered(:)) max(ref_filtered(:))])
ref_filt_plot = reshape(ref_filtered, 200,200,1,41);
montage(ref_filt_plot, 'DisplayRange', [min(ref_filt_plot(:)) max(ref_filt_plot(:))])
%%

figure(3)
montage(subt, 'DisplayRange', [min(subt(:)) max(subt(:))])





%% Classify


alpha = 0:dalpha:359;
n_rot = length(alpha); % number of rotations

best_ref = zeros(N_img,1);
best_rot = zeros(N_img,1);
for i=1:50 %N_img % loop through images
    disp(['Image ' num2str(i) ' of ' num2str(N_img)])
    % generate rotations of images and mirrors
    tic
    rotations = zeros(N, N, 2*n_rot, 'single');
    for j=1:n_rot
        rotations(:,:,j) = imrotate(img(:,:,i), alpha(j), 'crop');
    end
    for j=1:n_rot
        rotations(:,:,j+n_rot) = flip(imrotate(img(:,:,i), alpha(j), 'crop'),1);
    end
    toc
    tic
      %  figure(1)

    correlations = zeros(2*n_rot, N_ref);
    for r=1:2*n_rot % loop through rotations
        % calculate cc against each ref and retun maximum correlation +
        % index
       % disp(['Image ' num2str(i) ' of ' num2str(N_img) ', rotation ' num2str(r) ' of ' num2str(2*n_rot)])
        tmp = cc_against_classes(rotations(:,:,r), ref);
        correlations(r, :) = tmp'; % keep track of all correlations
        %    imagesc(correlations), axis image
        %    pause(.1)

    end
    toc
    [best_corr, i_tmp] = max(correlations(:));
    [best_rot(i), best_ref(i)] = ind2sub(size(correlations),i_tmp);

    
    
end

%%


for i=2:2 %N_img % loop through images

    figure(1)
    subplot(1,4,1)
    imagesc(img(:,:,i)), axis image
    title('Original image')
 
    subplot(1,4,2)
    mirrored = floor(best_rot(i)/n_rot);
    if mirrored
        rot_found =  flip(imrotate(img(:,:,i), alpha(best_rot(i)-n_rot), 'crop'),1);
    else
        rot_found = imrotate(img(:,:,i), alpha(best_rot(i)), 'crop');
    end
    imagesc(rot_found), axis image
    title('Rotation found')
    
    subplot(1,4,3)
    imagesc(ref_filtered(:,:,best_ref(i))), axis image
    title('Reference')

    x_cor = normxcorr2(ref_filtered(:,:,best_ref(i)), rotations(:,:,best_rot(i))) ;

    subplot(1,4,4)
    imagesc(x_cor(100:300,100:300) ), axis image, colorbar
    title('x-correlation image')
   pause
end


