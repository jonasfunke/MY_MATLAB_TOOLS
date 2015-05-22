function [ r_mean, sigma ] = dsDNA_mean_end_to_end_distances_semiflexible_polymer( L, exponent, lp)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    if nargin < 3
        lp = 150; % bp
    end
    dL = 0.001;

    r_mean = zeros(length(L),1);
    sigma = zeros(length(L),1);
    %r_mean3 = zeros(size(L));
    distribution = cell(length(L),1);

    for i=1:length(L)
        if L(i) == 0
            r_mean(i) = 0;
        else
            R = 0:dL:L(i)-dL;
            x = lp*(1-R./L(i))./L(i);
            %end_to_end = (4*pi*R.^2*(lp/L(i)^2).*frey_function(x)).^exponent;
            end_to_end = ((lp/L(i)^2).*frey_function(x)).^exponent;
            end_to_end = end_to_end./sum(end_to_end);
            distribution{i} = [R', end_to_end'./dL];
%             cumulant = zeros(length(R),1);
%             for j=1:length(cumulant)
%                cumulant(j) = sum(end_to_end(1:j));
%             end
% 
%             plot(R, end_to_end), hold on
%             tmp = find(cumulant-0.25 > 0);
%             vline(R(tmp(1)));
%             tmp = find(cumulant-0.75 > 0);
%             vline(R(tmp(1)));
          %   plot(R, end_to_end), hold on
          %   vline(sum(end_to_end.*R));
            r_mean(i, 1) = sum(end_to_end.*R); % mean
            sigma(i, 1) = sqrt( sum(end_to_end.*(r_mean(i, 1)-R).^2) ); % std
        end
    end
    
    

end

