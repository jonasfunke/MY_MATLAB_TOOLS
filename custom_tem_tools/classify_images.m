function [ classification ] = classify_images( img, references, dalpha )
% Sorts images img into references classes based on the maximum correlation
% between transformed (rotated and mirrored) image and reference class
%   Usage: out = classify_images(img, ref, 10);
    N_img = size(img,3);
    classification = zeros(N_img, 4);  % best_ref, alpha, mirror, cc
    parfor i=1:N_img % loop through images
        classification(i,:) = find_best_class(references, img(:,:,i), dalpha);
    end
    classification = [[1:N_img]' classification]; % add image index

end

