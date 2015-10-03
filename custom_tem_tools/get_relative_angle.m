function [ dalpha ] = get_relative_angle(varargin)
%Load data from imageJ(aquired with the arrow tool) and combutes relative
%angle
%%
if isempty(varargin)
    r_slice = 8;
    r_angle = 7;
else
    r_angle = varargin{1};
    r_slice = varargin{2};
end
%%
[fname, pname] =uigetfile('*.*', 'Select output file of imageJ.');
data = dlmread([pname fname], '\t', 1); % load data
%%
%keyboard
% convert absolute angles to relative angles

dalpha = zeros(size(data,1)/2,2);
for i=1:2:size(data,1)-1
    if data(i,r_slice) == data(i+1,r_slice)
        dalpha((i+1)/2,1) = abs(data(i+1,r_angle)-data(i,r_angle));
        dalpha((i+1)/2,2) = data(i,r_slice); 
        if dalpha((i+1)/2,1) > 180
            dalpha((i+1)/2,1) = 360-dalpha((i+1)/2,1);
        end
    else
        disp(['Warning: Out of slice sync at ' num2str(i) ])
    end
end
%%

if sum(diff(dalpha(:,2)==0)) > 0
    disp('Warning: some particles in wrong order')
    dalpha = sortrows(dalpha,2);
end

dlmwrite([pname fname(1:end-4) '_angles.txt'], dalpha, '\t') % write output

cf = figure();
xhist = 0:5:max(dalpha(:,1));
n = hist(dalpha(:,1), xhist);
bar(xhist, n)
pause
close(cf)
pause(0.1)

end

