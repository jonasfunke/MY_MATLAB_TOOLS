function [NN] = get_NN_density(xy,radius)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    x = xy(:,1);
    y = xy(:,2);
    NN = zeros(size(x));
    disp('computing density... please wait')
    for k = 1:length(x)
        v = (x-x(k)).^2 + (y-y(k)).^2 < radius^2;        % or sphere
        NN(k) = sum(v)-1;                           % number of points except x(i)
    end
    disp('density computed')

end

