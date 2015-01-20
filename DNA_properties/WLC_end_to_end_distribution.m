function [ p ] = WLC_end_to_end_distribution(t, r)
%Uprobability distribution of ssDNA from 
% Probing Single-Stranded DNA Conformational Flexibility Using Fluorescence
% Spectroscopy, M. C. Murphy,* Ivan Rasnik,* Wei Cheng,y Timothy M.
% Lohman,y and Taekjip Ha*, 2004
%   t = L / L_p; L = contour length, L_P = persistence length
%   r = R / L; R = end-to-end distance

% deterimin normalization
tmp = (3.*t./4);
A = 4 .* tmp.^(3/2) .* exp(tmp) ./ ( (4+12./tmp)+15./tmp.^2 ) ./  pi.^(3/2);

p = 4 .* pi .* A .* r.^2 .*  exp(-3.*t ./ (1-r.^2) ./ 4) ./ (1-r.^2).^(9/2);


end

