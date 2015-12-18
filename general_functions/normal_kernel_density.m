function [n, p, x_points] = normal_kernel_density(x_data, h, x_start, x_stop, dx)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    x_points = x_start:dx:x_stop;
    n = zeros(size(x_points));
    
    
    for j=1:length(x_data) % loop over data
        n = n + normpdf(x_points, x_data(j), h);
    end
    
    p = n/length(x_data);
%     x_hist = x_start:h:x_stop;
%     n_hist = hist(x_data, x_hist);
%     
%     cf = figure();
%     bar(x_hist, n_hist/sum(n_hist)/h), hold on
%     plot(x_points, p, 'r')
%     pause
%     close(cf)
end

