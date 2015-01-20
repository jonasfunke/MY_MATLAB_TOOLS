function [ R_mean ] = ssDNA_mean_end_to_end_distance( n_bases, lp)
% computes mean end-to-end distances based on 
% Probing Single-Stranded DNA Conformational Flexibility Using Fluorescence
% Spectroscopy, M. C. Murphy,* Ivan Rasnik,* Wei Cheng,y Timothy M.
% Lohman,y and Taekjip Ha*, 2004

    %lp = 2.3; % [bases] persitence length of ssDNA in units of bases


    R_mean = zeros(length(n_bases),1);
    for i=1:length(n_bases)
      %  t = n_bases(i) ./ lp;  % contour length / persitence length 
      %  p_norm = WLC_end_to_end_distribution(t, r) ; % propability
      %  end_to_end = p_norm ./sum(p_norm)/dr; % normalize
        if n_bases(i) == 0
             R_mean(i) = 0;
        else
              end_to_end = ssDNA_end_to_end_distribution( n_bases(i), lp );
              R_mean(i) = sum(end_to_end(:,2).*end_to_end(:,1))./sum(end_to_end(:,2));
        end

       % R = r.*n_bases(i); % actual end-to-end-distance in bases
       % R_mean(i) = sum(end_to_end.*R)./sum(end_to_end);     
    end
    
    

end

