function [n, p, x_points, ci] = uniform_kernel_density(x_data, h, x_start, x_stop, dx)
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
    p = n/length(x_data); % probability density
    
    % estimate errors using wilson score interval (this is the improved method) for binomial distribution confidence
    % interval 
    p_tmp =  p*h; % probability
    N = length(x_data); % number of total data points
     z = 1-0.05/2; % confidence interval 95%
    p_tmp_estimate = (p_tmp+z^2./2./N) ./ (1+z^2./N);
    p_tmp_err = (z./(1+z^2./N)).*sqrt( p_tmp.*(1-p_tmp)./N+z^2./4./N.^2 );
    ci = [p_tmp_estimate-p_tmp_err;  p_tmp_estimate+p_tmp_err]/h; % error in probabilty density
    
%     x_hist = x_start:h:x_stop;
%     n_hist = hist(x_data, x_hist);  
%     cf = figure();
%     bar(x_hist, n_hist/sum(n_hist)/h), hold on
%     plot(x_points, p, 'r')
%     pause
%     close(cf)
end

