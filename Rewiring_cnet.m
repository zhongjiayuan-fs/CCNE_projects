function [XY_edge] = Rewiring_cnet(x, y, k)
%DMI: compute the direct entropy causality from x to y.
% x: N samples
% y: N samples
% p: the number of neighbors of x 
% k: the k-th nearest number to use in calculating entropy, at least 2.
%Use the kNN method to estimate the MI(Y(xNN),YNN), where the dimension of YNN is p.
% define the neighbors of x 
% xi should be excluded from the neighbors of xi
N = length(x);
if nargin < 3
    k = 3; % k-th nearest neighbor
end
% define the neighbors of x 
% xi should be excluded from the neighbors of xi
dx = pdist2(x, x,'euclidean');
dy = pdist2(y,y,'euclidean');



[~,idx]=sort(dx,2);
[~,idy]=sort(dy,2);
XNNid= idx(:,1:k+1);
YNNid= idy(:,1:k+1);
XNNid(:,1)=1:N;
YNNid(:,1)=1:N;
YxNN = y(XNNid); % N*p 
YNN = y(YNNid);  % N*p

yfirst_column = YNN(:, 1);
YNN_dis=zeros(N,k);
YNN_dis(:, 1:end) = abs(YNN(:, 2:end) - yfirst_column);
[halfepsilon,halfepsilon_index]=max(YNN_dis,[],2);


xfirst_column = YxNN(:, 1);
YxNN_dis=zeros(N,k);
YxNN_dis(:, 1:end) = abs(YxNN(:, 2:end) - xfirst_column);

common_matrix = bsxfun(@lt,YxNN_dis,halfepsilon);
common_num=sum(common_matrix,2);

c=common_num/N;
a=(k * ones(N, 1))/N;
b=(k * ones(N, 1))/N;
XY_edge=(c.*(log(c./(a.*b))));
XY_edge(XY_edge<0)=0;
XY_edge = fillmissing(XY_edge,'constant',0);
end






