function [h ellipseXY] = plotEllipse(MU,COV,CI,dim,linewidth,linecolor,bLog)

% function [h ellipseXY] = plotEllipse(MU,COV,CI,dim,linewidth,linecolor,bLog)
%
%   example call: [h xy] = plotEllipse([0 0],[4 0; 0 16],95,[1 2],1,'k')
%
% plots ellipse for Gaussian w. MU and COV at specified confidence interval
%
% MU:        q-D mean vector       for nCtg levels [ nCtg x nF ]
% COV:       q-D covariance matrix for nCtg levels [ nF x nF x nCtg ]
% CI:        confidence internal
%            CI% of probability mass inside contour
% dim:       dimensions to plot                      [ 1 x 2 ]
% linewidth: duh!
% linecolor: duh!
% bLog:      whether MU and COV are in linear or log space
%            0 -> linear (default)
%            1 -> log
%%%%%%%%%%%%%%%%%%%%%
% h:         plot handle
% ellipseXY: x,y coordinates of ellipse

% INPUT HANDLING
if isempty(MU),                         MU = zeros(size(COV,3),size(COV,1));                                   end
if size(MU,2) ~= size(COV,1)            error(['plotEllipse: WARNING! dim of MU size must equal dim of COV']); end
if size(MU,1) == 1,                     MU = repmat(MU,size(COV,3));                                           end
if ~exist('dim','var') || isempty(dim)  dim = [1 2];                                                           end
if dim == 1,                            dim = [1 2];                                                           end
if ~exist('linewidth','var')            linewidth = 2;                                                         end
if ~exist('linecolor','var')            linecolor = repmat(get(gca,'colorOrder'),[5 1]);
                                        linecolor = linecolor(max([1 get(gca,'colororderIndex')-1]),:);        end
if ~exist('bLog','var')                 bLog = 0;                                                         end

%%%%%%%%%%%%%%%%%
% PLOT ELLIPSES %
%%%%%%%%%%%%%%%%%
for d = 1:size(COV,3)
    % EIGEN VECTORS AND VALUES OF COVARIANCE MATRIX
    [eigvec, eigval] = eig(COV(dim,dim,d));

    % LARGEST EIGEN VALUE
    eigvalMax       = max(diag(eigval));
    % LARGEST EIGEN VALUE INDEX
    [~,indCol] = find(eigval == eigvalMax);
    eigvecMax  = eigvec(:, indCol);

    % SMALLEST EIGEN VALUE
    eigvalMin  = min(diag(eigval));
    % SMALLEST EIGEN VALUE
    [~,indCol] = find(eigval == eigvalMin);
    eigvecMin  = eigvec(:, indCol);

    % ANGLE BETWEEN X-AXIS AND LARGEST EIGEN VECTOR
    rotDeg = mod(atan2d(eigvecMax(2), eigvecMax(1)),360);

    % ROTATION MATRIX
    R = [ cosd(rotDeg) sind(rotDeg); -sind(rotDeg) cosd(rotDeg) ];

    % ELLIPSE VALUE CORRESPONDING TO CONFIDENCE INTERVAL
    chi2val = sqrt(chi2inv(CI./100,2));

    % ELLIPSE PARAMETERS
    a = chi2val*sqrt(eigvalMax);
    b = chi2val*sqrt(eigvalMin);

    % PARAMETRIC EQUATION FOR ELLIPSE ( ellipse = a*cos(x) + b*sin(x)
    ellipseXY = bsxfun(@plus,MU(d,dim),[a*cosd( linspace(0,360,101)' ) b*sind( linspace(0,360,101)' )] * R);
    if bLog == 1
        ellipseXY = exp(ellipseXY);
    end
    % PLOT ELLIPSE
    h = plot(ellipseXY(:,1),ellipseXY(:,2),'linewidth',linewidth,'color',linecolor);
    hold on
end
