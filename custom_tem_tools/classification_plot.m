%%

% load file with data

% load classification file

% 
% for i=100:200%N_img % loop through images
% 
%     figure(1)
%     subplot(1,4,1)
%     imagesc(img(:,:,i)), axis image
%     title(['Original image: ' num2str(i)]  )
%  
%     subplot(1,4,2)
%     rot_found = my_imagetransform(img(:,:,i), classification(i,3), classification(i,2)); % mirror, angle
%     imagesc(rot_found), axis image
%     title('Rotated image')
%     
%     subplot(1,4,3)
%     imagesc(ref_filtered(:,:,classification(i,1))), axis image
%     title(['Reference ' num2str(classification(i,1))])
% 
%     x_cor = normxcorr2(ref_filtered(:,:,classification(i,1)), rot_found) ;
%     x_cor = x_cor(100:300,100:300);
%     subplot(1,4,4)
%     hold off
%     imagesc(x_cor ), axis image, hold on
%     [tmp_cor, i_tmp] = max(x_cor(:)); % use best overall correlation
%     [best_i, best_j] = ind2sub(size(x_cor),i_tmp);
%     plot(best_j, best_i, 'ro')
%     title(['x-correlation image (cor = ' num2str(tmp_cor) ')'])
%     
%     figure(2)
%     hold off
%     imagesc(cor_matrix{i}), axis image, colorbar, hold on
%  
%     [~, i_tmp] = max(cor_matrix{i}(:)); % use best overall correlation
%     [best_rot, best_ref] = ind2sub(size(cor_matrix{i}),i_tmp);
%     plot(best_ref, best_rot, 'ro')
%    
%     figure(3)
%     subplot(1,2,1)
%     hold off
%     tmp =min(cor_matrix{i});
%     bg = tmp(ones(144,1),:);
%     imagesc(cor_matrix{i}-bg), axis image, colorbar, hold on
%     
%     [~, i_tmp] = max(cor_matrix{i}(:)-bg(:)); % use best overall correlation
%     [best_rot, best_ref] = ind2sub(size(cor_matrix{i}),i_tmp);
%     plot(best_ref, best_rot, 'ro')
%     
%     subplot(1,2,2)
%     imagesc(ref_filtered(:,:,best_ref)) , axis image
%     
%     
%    pause
% end
% 

%% correct classification for background

alpha = 0:dalpha:359;  % angles
n_rot = length(alpha); % number of rotations
classification_cor = zeros(size(classification,1), 4);
cor_matrix_cor = cell(size(cor_matrix));

for i=1:1000 % loop through images

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
    classification_cor(i,:) = [best_ref alpha_found mirror max_cor]; % combine output
    
end 



for i=1:200%N_img % loop through images
    figure(1)
    
    subplot(2,2,1)
    rot_found = my_imagetransform(img(:,:,i), classification_cor(i,3), classification_cor(i,2)); % mirror, angle
    imagesc(rot_found), axis image
    title(['Rotated image ' num2str(i)])
    
    subplot(2,2,2)
    imagesc(ref_filtered(:,:,classification_cor(i,1))), axis image
    title(['Reference ' num2str(classification_cor(i,1)) ', (cor = ' num2str(tmp_cor) ')'])

   subplot(2,2,[3 4])

    hold off
    imagesc(cor_matrix_cor{i}'), axis image, colorbar, hold on
    
    [~, i_tmp] = max(cor_matrix_cor{i}(:)); % use best overall correlation
    [best_rot, best_ref] = ind2sub(size(cor_matrix_cor{i}),i_tmp);
    plot( best_rot,best_ref, 'ro')
   
    xlabel('Rotation'), ylabel('Reference')
    title('Maximum correlation')
   pause
end





