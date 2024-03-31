classdef AMANrn < handle & AMAChld
%ctgInd  [nStm x 1]
%X       [1    x nUCtg]
%s       [nPix  x nStm]     o
%f       [nPix  x nF]       o
%r       [nStim x nF]]      x
properties
    f
    w
    r
    R
    Rm

    bFourierF % ??? MOVE MAYBE?
    fCmplx
    bAnalytic
    % cmplx as
    %    0 single response
    %    1 single response XXX
    %    2 two responsese (real & imag)
    %    3 two responsese (abs & arg)

    rMax
    RMaxType
    fano

    % noise
    bNoise
    var0
    rho
    normType

    nF
    nvar
    ns
    Lambda
end
properties(Hidden)
    stimType
    bUpdate=true

    AMADeps={'Stim'}
end
properties(Access=private)
    Stim
end
methods(Static)
    function P=getP()
        constTest= @(x) isempty(x) | (numel(x)==1 & isnumeric(x));
        P={'rMax',5.7,constTest;
           'fano',1.36,constTest;
           'var0',0.23,constTest;
           'rho',0,constTest;
           ...
           'bFourierF',false,'isBinary_e';
           'fCmplx',[],'isBinary_e';
           'bAnaltytic',[],'isBinary_e';
           ...
           'bNoise',false,'isBinary_e';
           'normType','none','ischar_e';
           'RMaxType',0,'Num.is';
       };
    end
end
methods
    function obj=AMANrn()
        if nargin > 1
            return
        end
    end
    function init(obj)
        if isempty(obj.fCmplx)
            if ~isempty(obj.bAnalytic)  && obj.bAnalytic > 0
                obj.fCmplx=1;
            else
                obj.fCmplx=double(obj.bFourier);
            end
        end
        if isempty(obj.bAnalytic)
            obj.bAnalytic = obj.fCmplx > 0;
        end
        obj.update();
    end
    function assert(obj)
        assert(isempty(obj.f0) || size(obj.f0,1) == obj.nPix,...
               'size of f0 dim 1 (%d) does not match number of stim (%d).',size(obj.f0,1),obj.nPix);
    end
    function out=get.Stim(obj)
        out=obj.AMA.Stim;
    end
%- SET
    function update(obj)
        if ~obj.bUpdate
            return
        end

        % INIT W/ STIM
        S=obj.Stim;
        if obj.bFourierF && strcmp(classType,'AMAStim')
            S.bFourierS=true;
            if isempty(S.As)
                S.getAmplitudes();
            end
            % XXX MAKE SURE S DIMENSIONS ARE CORRECT
            % obj.S.reshape
        end
    end
    function set.f(obj,f)
        obj.nF=size(f,2);
        obj.f=f;
    end
    function set.Stim(obj,S)
        obj.stimType=class(S);
        obj.Stim=S;
        obj.bUpdate=true;
    end
%- GET
    function respond(obj,S)
        obj.update();
        if nargin < 2 || isempty(S)
            S=obj.Stim;
        else
            obj.stimType=class(S);
        end
        % RESPOND
        obj.r =(S*obj.f);

        % NORMALIZE
        Rm=obj.normalize(S);

        % NOISE
        obj.getNoise();
        if obj.bNoise
            nu=obj.ns;
        else
            nu=0;
        end

        obj.R  = obj.rMax * (Rm+nu);
        obj.Rm = obj.rMax *  Rm;
    end
    function R=normalize(obj,S)
        switch obj.normType
        case 'nrw'
            if strcmp(obj.stimType,'AMAStim')
                R=obj.r./(S.As' * abs(obj.f));
                %N=dot(abs(obj.f),abs(S.As));
            else
                error('not implemented') % TODO
            end
        case 'brd'
            if strcmp(obj.stimType,'AMAStim')
                R=obj.r./sqrt(sum(S.S.^2,2));
            end
        case {'none',''}
            R=obj.r;
            return
        otherwise
            error('Invalid norm type')
        end

        % MAXR
        switch obj.fCmplx
        case 0
            switch obj.RMaxType
            case 1
                R=R./max(R(:))*max(obj.r(:));
            case 2
                R=R./max(R(:));
            case 3
                ind=abs(R)>1;
                R(ind)=1*sign(R(ind));
            end
        case 1
            switch obj.RMaxType
            case 1
                % rescale to max r
                if obj.fCmplx
                R=R./obj.cmpMaxAbs(R(:),2)*obj.cmpMaxAbs(obj.r(:),2);
            case 2
                % rescale to 1
                R=R./obj.cmpMaxAbs(R(:),2);
            case 3
                % truncate over 1
                ind=obj.cmpMaxAbs(R(:),2)>1;
                R(ind)=1*cmpSgn(R(ind));
            end
        case 2
            % XXX
        end
    end
    function getNoise(obj)
        switch obj.fCmplx
        case 0
            obj.nvar=sqrt( obj.fano .* abs(obj.R) + obj.var0 );
        case 1
            obj.nvar=sqrt( obj.fano .* obj.cmpAbs(obj.R(:)) + obj.var0 );
        case 2
            % XXX
        end
        obj.ns = normrnd(0,obj.nvar,size(obj.R));
    end
    function getLambda(obj)
        % noise variance
        % MATCH CONSTANT ADDITIVE TO AVERAGE SCALED ADDITIVE NOISE VARIANCE
        varAvg     = mean(obj.nvar.^2);
        % INTERNAL FILTER RESPONSE COVARIANCE MATRIX (ASSUMING UNCORRELATED NOISE)
        obj.Lambda     = diag( varAvg.*ones(1,obj.nF) );

        if obj.rho==0
            return
        end

        % CORRELATED INTERNAL NEURAL NOISE %%
        C = corr2cov(varAvg,obj.rho);
        obj.Lambda(logical(triu(ones(obj.nF), 1))) = C;
        obj.Lambda(logical(tril(ones(obj.nF),-1))) = C;

    end
%- UTIL
    function out=disp_maxR(obj)
        str=sprintf('  Max R       = %.3f',max(obj.R(:)));
        if nargout > 0
            out=str;
        else
            disp(str);
        end
    end
end
methods(Static)
    function out=cmpSgn(R)
        out=sign(real(R)) + i*signIm(R);
    end
    function out=cmpMaxAbs(R,dim)
        if nargin < 2 || isempty(dim)
            dim=3;
        end
        max( abs( cat(dim,real(R),imag(R)) ) ,[],dim);
    end
    function out=cmpAbs(R,dim)
        if nargin < 2 || isempty(dim)
            dim=3;
        end
        out=abs(mean( cat(dim,real(R),imag(R)) ,dim));
    end
    function f0=initRandomFilters(f0,nPix,nF,rndSd)
        if nargin < 4 || isempty(rndSd)
            rndSd=1;
        end
        rng('default');

        rng(rndSd);
        %nStim=size(s,1);
        % RANDOM FILTER WEIGHTS
        f0 = [f0(:,1:min([nF size(f0,2)])) randn(nPix,nF-size(f0,2))];
        % NORMALIZED TO L2 NORM OF 1.0
        f0 = bsxfun(@rdivide,f0,sqrt(sum(f0.^2,1)));
    end
    function w0=initRandomWeights(w0,nF,rndSd)
        if nargin < 3 || isempty(rndSd)
            rndSd=1;
        end
        rng('default');
        rng(rndSd);
        %nStim=size(s,1);
        % RANDOM FILTER WEIGHTS
        %w0 = [w0(:,1:min([nF size(w0,2)])) rand(nF,nF-size(w0,2))];
        w0 = [w0(:,1:min([nF size(w0,2)])) ones(nF,nF-size(w0,2))];
        % NORMALIZED TO L2 NORM OF 1.0
        w0 = bsxfun(@rdivide,w0,sqrt(sum(w0.^2,1)));
    end
end
end
