function [ r_mean ] = dsDNA_mean_end_to_end_distances_semiflexible_polymer( L, exponent )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    lp = 50; % bp
    dL = 0.001;

    r_mean = zeros(size(L));
    %r_mean3 = zeros(size(L));

    for i=1:length(L)
        if L(i) == 0
            r_mean(i) = 0;
        else
            R = 0:dL:L(i)-dL;
            x = lp*(1-R./L(i))./L(i);
            end_to_end = (4*pi*R.^2*(lp/L(i)^2).*frey_function(x)).^exponent;
            r_mean(i) = sum(end_to_end.*R)./sum(end_to_end);
        end
    end

end

