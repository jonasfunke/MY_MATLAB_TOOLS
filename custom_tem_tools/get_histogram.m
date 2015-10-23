function [ alpha, n, p ] = get_histogram( data, dalpha )
% computes the histogram of of a data struct 

    alpha = 0:dalpha:120;
    n = hist(data.angles, alpha);
    n = n';
    p = n./sum(n)./dalpha;

end

