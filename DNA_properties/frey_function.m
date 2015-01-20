function [ y ] = frey_function( x )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
y = zeros(size(x));
for i=1:length(x)
    if x(i) > 0.2
        y(i) = pi.*exp(-pi.^2.*x(i))/2;
    else
        y(i) = (1./x(i)-2) .* exp(-1./x(i)/4)./(8.*pi^1.5.*x(i).^1.5);
    end
end

end

