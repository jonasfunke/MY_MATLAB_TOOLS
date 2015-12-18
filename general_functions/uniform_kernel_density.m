function [n, p, x_points] = uniform_kernel_density(x_data, h, x_start, x_stop, dx)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    x_points = x_start:dx:x_stop;
    n = zeros(size(x_points));
    for i=1:length(x_points)
        x_cur = x_points(i);
        n(i) = (0.5/h).*sum((x_data > x_cur-h) & (x_data < x_cur+h));
%         for j=1:length(x_data)
%             if (x_data(j)-h/2 < x_points(i)) && (x_points(i) <= x_data(j)-h/2)
%                 n(i) = n(i)+1;
%             end
%         end
    end
    p = n/length(x_data);
    
    
%     x_hist = x_start:h:x_stop;
%     n_hist = hist(x_data, x_hist);  
%     cf = figure();
%     bar(x_hist, n_hist/sum(n_hist)/h), hold on
%     plot(x_points, p, 'r')
%     pause
%     close(cf)
end

