function colors = getColorOrder(co,numColors)

% function colors = getColorOrder(co,numColors)
%
%   example call: colors = getColorOrder([],length(X))
%
% co:        color order (or color map)
% numColors: number of colors in the colormap
% %%%%%%%%%%%%%
% colors:    color order of length numColors


if ~exist('co','var') || isempty(co)
    co  = [ 0         0    1.0000;
            .25    0.5000        0.25;
            1.0000         0         0;
            0           0.7500    0.7500;
            0.7500         0    0.7500;
            0.7500    0.7500         0;
            0.2500    0.7500    0.2500];
    co      = [flipud(co(2:end,:)); co];
end

colors  = interp1([1:size(co,1)]',co,linspace(1,size(co,1),numColors));