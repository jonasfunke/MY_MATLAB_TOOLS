function [ output_args ] = find_best_class( ref, img )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


    disp(['Image ' num2str(i) ' of ' num2str(N_img)])
    % generate rotations of images and mirrors
    
    rotations = zeros(N, N, 2*n_rot, 'single');
    for j=1:n_rot
        rotations(:,:,j) = imrotate(img(:,:,i), alpha(j), 'crop');
    end
    for j=1:n_rot
        rotations(:,:,j+n_rot) = flip(imrotate(img(:,:,i), alpha(j), 'crop'),1);
    end
    
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
    
   % [best_corr, i_tmp] = max(correlations(:));
   % [best_rot(i), best_ref(i)] = ind2sub(size(correlations),i_tmp);
end

