function label = legendLabel(textPrefix,numericalPostFix,bExpandScalar,strLength)

% function legendLabel(textPrefix,numberOfLabels)
% 
%   example call: % SIMPLE LABEL IN UPPER RIGHT (i.e. 'northeast') WITH NO DECIMALS 
%                 legend( legendLabel('n=',n,1,0),'location','northeast') );
%
%                 % SIMPLE LABEL IN LOWER RIGHT (i.e. 'southeast') WITH TWO DECIMALS 
%                 legend(legendLabel('n=',n,1,2),'location','southeast'));
%
% spits out cell array of strings to be used with legend.m
% 
% textPrefix:     string containing text to precede number
% numberOfLabels: a scaler or 1xn vector
% bExpandScalar:  
% strLength:      numstr('string',strlength). default = 3

if ~exist('bExpandScalar','var') || isempty(bExpandScalar)
    bExpandScalar = 1;
end
if ~exist('strLength','var') || isempty(strLength)
    strLength = 2;
end

if length(numericalPostFix) > 1 || bExpandScalar == 0
    label = [repmat(textPrefix,length(numericalPostFix),1) num2str([numericalPostFix]',['%.' num2str(strLength) 'f'])];
else
    label = [repmat(textPrefix,length(numericalPostFix),1) num2str(numericalPostFix',['%.' num2str(strLength) 'f'])];
end

label = cellstr(label);
