function d =  circ_dist(x,y)
%
% d = circ_dist(alpha, beta)
%
%   Pairwise difference x_i-y_i around the circle computed efficiently.
%
%   Input:
%     x       sample 1 of circular variable
%     y       sample 2 of circular variable
%
%   Output:
%     d       matrix with distances
%
% References:
%     Biostatistical Analysis, J. H. Zar, p. 651
%
% PHB 3/19/2009
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html


if size(x,1)~=size(y,1) && size(x,2)~=size(y,2) && length(y)~=1
  error('Input dimensions do not match.')
end

% EASY TO READ (LESS FLEXIBLE)
% d = angle(exp(1i*x)./exp(1i*y));

% HARDER TO READ (MORE FLEXIBLE)
d = angle(bsxfun(@rdivide,exp(1i*x),exp(1i*y)));