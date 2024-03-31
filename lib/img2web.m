function [Iweb,DC] = img2web(I,W,bPreWndw,nChnl)

% function [Iweb,DC] = img2web(I,W,bPreWndw,nChnl)
%
%   example call: img2web(I,W,bPreWndw,nChnl)
%
% convert intensity image to weber contrast image
%
% I:          intensity image                              [ Nd x nStm ]
% W:          window                                       [ Nd x  1   ]
% bPreWndw:   boolean indicating if image is pre-windowed
%             1 -> prewindowed
%             0 -> not so much
% nChnl:      number of channels in which the normalization should take place
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iweb:       weber contrast image                         [ Nd   x nStm  ]
% DC:         DC                                           [ nStm x nChnl ]

tol = 1e-1;

% INPUT HANDLING
if sum(mean(I) < tol)>1 error(['img2web: WARNING! ' num2str(sum(mean(I)<tol)) ' of the I images are near mean = 0. I may be a contrast image, though it should be an intensity image']); end
if ~exist('W','var')        || isempty(W)           W = ones(size(I,1),1); disp(['img2web: WARNING! no windowing will be performed']);   end
if ~exist('bPreWndw','var') || isempty(bPreWndw)    error(['img2web: WARNING! bPreWndw must be defined (0 or 1) to use img2web.m ...']); end
if bPreWndw ~= 0 && bPreWndw ~= 1                   error(['img2web: WARNING! bPreWndw must be defined (0 or 1) to use img2web.m ...']); end
if ~exist('nChnl','var')    || isempty(nChnl)       nChnl   = 1;                                                                         end

% CHANNEL INDICES TO NORMALIZE TOGETHER [ Nd/nChnl x nChnl ]
indNrm   = reshape(1:size(I,1),[],nChnl);

% PREALLOC MEMORY
Iweb = zeros(size(I));
if bytes(I) <= 1e9
    for c = 1:nChnl
        % WEBER CONTRAST IMAGE FOR EACH CHANNEL
        [Iweb(indNrm(:,c),:),DC(c,:)] = contrastImageVec(I(indNrm(:,c),:),W(indNrm(:,c)),bPreWndw);
    end
else
   for c = 1:nChnl 
        for i = 1:size(I,2)
        [Iweb(indNrm(:,c),i),DC(c,i)] = contrastImage(I(indNrm(:,c),i),W(indNrm(:,c)),bPreWndw);
        end
   end
end
DC = DC';

