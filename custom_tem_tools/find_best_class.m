function [ output ] = find_best_class(ref, img, dalpha)
% determines in which class (of references ref) the image img belongs
%   Detailed explanation goes here

    alpha = 0:dalpha:359;  % angles
    n_rot = length(alpha); % number of rotations

    % generate rotations of images and mirrors
    rotations = zeros(size(img,1), size(img,1), 2*n_rot, 'single'); % include mirror transformations
    for j=1:n_rot
        rotations(:,:,j) = my_imagetransform(img, 0, alpha(j)); % no mirror
    end
    for j=1:n_rot
        rotations(:,:,j+n_rot) = my_imagetransform(img, 1, alpha(j)); % mirror
    end
    
    correlations = zeros(2*n_rot, size(ref,3)); % matrix that stores max correlation with each class
    for r=1:2*n_rot % loop through rotations
        % calculate cc against each ref and retun maximum correlation
        tmp = cc_against_classes(rotations(:,:,r), ref);
        correlations(r, :) = tmp'; % keep track of all correlations
    end
    
    [best_corr, i_tmp] = max(correlations(:)); % use best overall correlation
    [best_rot, best_ref] = ind2sub(size(correlations),i_tmp);
    
    disp(['Found reference: ' num2str(best_ref) ', Found rotation: ' num2str(best_rot) ])
    mirror = best_rot > n_rot; % determine if image was mirrored
    % determine rotation angle
    if mirror 
        alpha_found = alpha(best_rot-n_rot);
    else
        alpha_found = alpha(best_rot);
    end
    output = [best_ref alpha_found mirror best_corr]; % combine output
    
    
    % plot found image
%     subplot(1,3,1)
%     imagesc(img), axis image
%     title('Image')
%     
%     subplot(1,3,2)
%     imagesc(rotations(:,:,best_rot)), axis image
%     title('Rotated image')
%     
%     subplot(1,3,3)
%     imagesc(ref(:,:,best_ref)), axis image
%     title('Reference')
    
    
end

