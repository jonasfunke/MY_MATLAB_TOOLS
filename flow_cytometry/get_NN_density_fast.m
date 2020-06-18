function [NN, NN_tmp, xy] = get_NN_density_fast(xy,resolution, varargin)
% calculate local density based on a grid
%   Detailed explanation goes here

% remove data with NaN or Inf
i_discard = isnan(xy(:,1)) | isnan(xy(:,2)) | isinf(xy(:,1)) | isinf(xy(:,2));
i_use = ~i_discard;

if ~isempty(varargin)    
    limits = varargin{1};
    i_use2 = xy(:,1)>limits(1) & xy(:,1)<limits(2) & xy(:,2)>limits(1) & xy(:,2)<limits(2) ;
    i_use = i_use & i_use2;
end

xy = xy(i_use,:);

% 
% xy(isnan(xy(:,1)) | isnan(xy(:,2)),:) = [];
% xy(isinf(xy(:,1)) | isinf(xy(:,2)),:) = [];



Nx = ceil((max(xy(:,1))-min(xy(:,1)))/resolution);
Ny = ceil((max(xy(:,2))-min(xy(:,2)))/resolution);
if Nx==0
    Nx=1;
end
if Ny==0
    Ny=1;
end


% put data into grid
i = max(ceil( (xy(:,1)-min(xy(:,1)))/resolution),1); % this sets the min to 1
j = max(ceil( (xy(:,2)-min(xy(:,2)))/resolution),1);
xygrid = zeros(Nx, Ny, 'int16');
for k=1:size(xy,1)
    xygrid(i(k), j(k)) = xygrid(i(k), j(k)) +1;
end

% get the NN 
NN_tmp = xygrid(sub2ind(size(xygrid), i,j));

NN = zeros(length(i_discard), 1);
NN(i_discard) = 0;
NN(i_use) = NN_tmp;
end

