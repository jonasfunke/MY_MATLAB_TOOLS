function [ p, A_sub, B_sub] = calculateRation( A, B_in, status )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


    
    
    %[cc, shift, lane_da_shift] = xcorr2_bounded(lane_dd, lane_da, 5); % find best overlay of images
    [cc, shift, B] = xcorr2_bounded(A, B_in, 5, 0); % find best overlay of images
    
    % only use subimage for further analysis
    dy = shift(2);
    dx = shift(1);
    
    B_sub = B( max(1,1+dy):min(size(B,1), size(B,1)+dy), max(1,1+dx):min(size(B,2), size(B,2)+dx) );
    A_sub = A( max(1,1+dy):min(size(B,1), size(B,1)+dy), max(1,1+dx):min(size(B,2), size(B,2)+dx) );

    %
    %{
    %first algorith
    profileA = sum(A,2);
    profileB = sum(B,2);
    R(1,1) = max(profileA) ./ max(profileB);
    
    % second algorithm
    R(2,1) = sum(sum(A)) ./ sum(sum(B));
    
    %}
    
    %third algorithm
    p = polyfit(B_sub(:), A_sub(:), 1);
    p_raw = polyfit(B_in(:), A(:), 1);
    
    R = p(1);
    
    x = [min([B_sub(:); B_in(:)]) max([B_sub(:); B_in(:)])];
    
    
   % subplot(2,1,1)
   % plot(1:size(A,1), profileA, 'r', 1:size(B,1), profileB, 'g')
   % legend({'A', 'B'})
     
   if status
        %close all
        %{
        subplot(2,1,1)
        plot(B_in(:), A(:), 'b.', x, p_raw(1).*x+p_raw(2), 'r', x,  (sum(sum(A)) ./ sum(sum(B))).*x, 'g')
        legend({'raw data', ['raw fit: ' num2str(p_raw(1))], num2str(sum(sum(A)) ./ sum(sum(B))) })
        
        subplot(2,1,2)
        plot(B(:), A(:), 'b.', x, p(1).*x+p(2), 'r', x,  (sum(sum(A)) ./ sum(sum(B))).*x, 'g')
        legend({'shifted data', ['shifted fit: ' num2str(p(1))], num2str(sum(sum(A)) ./ sum(sum(B)))})
        %}
        
        plot(B_in(:), A(:), 'b.', x, p_raw(1)*x+p_raw(2), 'b', B_sub(:), A_sub(:), 'r.', x, p(1)*x+p(2), 'r')
        
   end
    
    
    

end

