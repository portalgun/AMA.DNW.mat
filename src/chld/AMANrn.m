classdef AMANrn < handle & AMAChld
%ctgInd  [nStm x 1]
%X       [1    x nUCtg]
%s       [nPix  x nStm]     o
%f       [nPix  x nF]       o
%r       [nStim x nF]]      x
properties
    f
    w

    r %  raw
    rc % raw combined
    Rn % normalized
    Rc % normalized, combined
    Rm % normalized, combined, attenuated
    R  % normalized, combined, attenuated, noised

    bFourierF % ??? MOVE MAYBE?

    bSplitR

    fCmplx
    % cmplx as
    %    0 single response (real)
    %    1 single response (real and imag)
    %    2 two responsese (real & imag)
    %    3 two responsese (real & real)

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
        P=[...
            AMANrn.getP_large();
            AMANrn.getP_base();
        ];
    end
    function P=getP_large()
        P={'f',[],'';
           'w',[],'';
        };
    end
    function P=getP_base()
        constTest= @(x) isempty(x) | (numel(x)==1 & isnumeric(x));
        P={'rMax',5.7,constTest;
           'fano',1.36,constTest;
           'var0',0.23,constTest;
           'rho',0,constTest;
           %coefs %*% xmat matrix multiplication
           ...
           'bFourierF',false,'isBinary_e';
           'bSplitR',false,'isBinary_e';
           'fCmplx',[],'';
           ...
           'bNoise',false,'isBinary_e';
           'normType','none','ischar_e';
           'RMaxType',0,'Num.is';
       };
    end
end
methods
    function P=getPSave(obj)
        P=[...
            obj.getP_base();
        ];
    end
    function obj=AMANrn()
        if nargin > 1
            return
        end
    end
    function init(obj)
        if isempty(obj.fCmplx)
            obj.fCmplx=double(obj.bFourierF);
        end

        obj.bUpdate=true;
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
        if isempty(obj.stimType)
            obj.stimType=class(obj.Stim);
        end
        if obj.bFourierF && strcmp(obj.stimType,'AMAStim')
            S.bFourierS=true;
            if isempty(S.As)
                S.getAmplitudes();
            end
            % XXX MAKE SURE S DIMENSIONS ARE CORRECT
            % obj.S.reshape
        end
        obj.bUpdate=false;
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
    function append(obj)
        obj.prps=[obj.prps; 'bSplitR'];
    end
    function respond(obj,S)
        if isempty(obj.stimType)
            obj.bUpdate=true;
        end
        obj.update();
        if nargin < 2 || isempty(S)
            S=obj.Stim;
        else
            obj.stimType=class(S);
        end
        % RESPOND
        if obj.fCmplx==1
            f=real(obj.f);
        else
            f=obj.f;
        end
        obj.r=mtimes(S,f,obj.bSplitR); % XXX SLOW 2

        % NORMALIZE
        obj.Rn=obj.normalize(obj.r,S); % XXX SLOW 3

        % Merge quadrature responses
        obj.rc=obj.combine(obj.r);
        obj.Rc=obj.combine(obj.Rn);

        % apply maximum FR
        obj.Rm=obj.attenuate(obj.Rc,obj.rc);

        % NOISE
        obj.getNoise(obj.Rm); % XXX SLOW 1
        if obj.bNoise
            nu=obj.ns;
        else
            nu=0;
        end

        obj.R  = obj.rMax * (obj.Rm+nu);
        obj.Rm = obj.rMax *  obj.Rm;
    end
    function Rn=normalize(obj,r,S)
        switch obj.normType
        case 'nrw'
            if strcmp(obj.stimType,'AMAStim')
                Rn=r./(S.As' * abs(obj.f));
                %N=dot(abs(obj.f),abs(S.As));
            else
                error('not implemented') % TODO
            end
        case 'brd'
            if strcmp(obj.stimType,'AMAStim')
                Rn=r./sqrt(sum(S.S.^2,2));
            end
        case {'none',''}
            Rn=r;
        otherwise
            error('Invalid norm type')
        end
    end
    function Rc=combine(obj,R)
        if obj.fCmplx==3
            Rc=[real(R) imag(R)];
        else
            Rc=R;
        end
    end
    function Rm=attenuate(obj,R,r)
        % MAXR
        switch obj.fCmplx
        case {0,1,3}
            switch obj.RMaxType
            case 0
                Rm=R;
            case 1
                Rm=R./max(R(:))*max(r(:));
            case 2
                Rm=R./max(R(:));
            case 3
                ind=abs(R)>1;
                Rm=R;
                Rm(ind)=1*sign(R(ind));
            end
        case 2
            switch obj.RMaxType
            case 0
                Rm=R;
            case 1
                % rescale to max r
                Rm=R./obj.cmpMaxAbs(R(:),2)*obj.cmpMaxAbs(obj.r(:),2);
            case 2
                % rescale to 1
                Rm=R./obj.cmpMaxAbs(R(:),2);
            case 3
                % truncate over 1
                ind=obj.cmpMaxAbs(R(:),2)>1;
                Rm=R;
                Rm(ind)=1*cmpSgn(R(ind));
            end
        end
    end
    function getNoise(obj,R)
        switch obj.fCmplx
        case {2,3}
            obj.nvar=sqrt( obj.fano .* obj.cmpAbs(R) + obj.var0 );
            % cmpAbs XXX SLOW 2
        case 3
            % XXX
            %
        otherwise
            obj.nvar=sqrt( obj.fano .* abs(R) + obj.var0 );
        end
        %obj.ns = obj.normrnd(obj.nvar,size(R));
        obj.ns = randn(size(R),'like',obj.nvar) .* obj.nvar;
    end

    function getLambda(obj)
        % noise variance
        % MATCH CONSTANT ADDITIVE TO AVERAGE SCALED ADDITIVE NOISE VARIANCE
        varAvg     = mean(obj.nvar.^2);
        % INTERNAL FILTER RESPONSE COVARIANCE MATRIX (ASSUMING UNCORRELATED NOISE)
        obj.Lambda     = diag( varAvg.*ones(1,size(obj.R,2)) );

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
    function out=splitF(obj,splitI,fI)
        % [nPix/nSplit nSplit nF]

        Stim=obj.AMA.Stim;
        if nargin < 2 || isempty(splitI)
            splitI=Stim.getSplitInds(1:Stim.nSplit);
        end
        if nargin < 3 || isempty(fI)
            fI=1:obj.nF;
        end
        nF=numel(fI);
        nPix=size(splitI,1);
        nI=size(splitI,2);

        f=obj.f(:,fI);
        out=zeros(nPix,nI,nF);
        for i = 1:nI;
            inds=splitI(:,i);
            for j = fI
                out(:,i,j)=f(inds,j);
            end
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
