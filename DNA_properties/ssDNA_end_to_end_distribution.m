function [ end_to_end_dist ] = ssDNA_end_to_end_distribution( n_bases, lp )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
   
 %  lp = 2.3; % [bases] persitence length of ssDNA
    
   t = n_bases ./ lp; % contour length / persitence length 
   dr = 0.001;
   r = [0:dr:1-dr];
   

   % if strcmp(mode, 'frey')
   %     x = (1-r)./t;
   %     end_to_end = 4*pi * r.^2 .* frey_function(x); % left lp out - only normalization
   %     end_to_end = end_to_end ./ sum(end_to_end)/dr;
   % else % use the TJ Ha model
        p_norm = WLC_end_to_end_distribution(t, r) ; % propability
        end_to_end = p_norm ./sum(p_norm)/dr; % normalize
   % end
    
    
    end_to_end_dist = [r'.*n_bases, end_to_end'./n_bases];
    
    
    
       
        
end

