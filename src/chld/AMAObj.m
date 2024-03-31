classdef AMAObj < handle & AMAChld & AMAmvn
%fMM X
%FLL X
%GSS X
%GRD
%SGD
properties
    fun

% OUTPUT
    LAll % [ nStim x nCtg ] @ correct category
    LCor % [ nStim x 1 ]

    PAll
    PCor

    XHat
    XDev

    E
    Em
    Em0

end
% XXX MAKE READ ONLY
properties
    f
    w

    Opts
    N
    S
    Index % filterIndex
    FMinOpts
end
properties(Hidden)
    AMADeps={'N','S','Index','FMinOpts'}
    bErrored

    ctgInd
    nStim
    nCtg
    nF
    mvncnst
    bRan=false
    lcorind

    aggIdx
    yagg0
end
methods
%- MAIN
    function obj=AMAObj()
        if nargin < 1
            return
        end
    end
    function out=get.Opts(obj)
        out=obj.AMA.Opts.Obj;
    end
    function out=get.N(obj)
        out=obj.AMA.Nrn;
    end
    function out=get.f(obj)
        out=obj.AMA.Nrn.f;
    end
    function out=set.f(obj,val)
        obj.AMA.Nrn.f=val;
    end
    function out=set.w(obj,val)
        obj.AMA.Nrn.w=val;
    end
    function out=get.S(obj)
        out=obj.AMA.Stim;
    end
    function out=get.Index(obj)
        out=obj.AMA.Index;
    end
    function out=get.FMinOpts(obj)
        out=obj.AMA.Opts.FMinOpts;
    end
%- UTIL
    function [n,N,c,bCorrect]=get_nCorrect(obj,errType)
        if nargin < 2 || isempty(errType)
            errType=obj.Opts.errType;
        end
        switch errType
        case 'ME'
            ex=1;
        otherwise
            ex=2;
        end
        err=obj.E;
        tol=abs((obj.S.X(2)-obj.S.X(1))/2);
        bCorrect=err < tol.^(1./ex);
        n=sum(bCorrect);
        N=length(err);
        c=n/N;
    end
    function out=disp_nCorrect(obj,errType)
        if nargin < 2
            errType=[];
        end
        [n,N]=obj.get_nCorrect(errType);
        str=sprintf('  Correct %5d / %d = %.3f',n,N,n/N*100);
        if nargout > 0
            out=str;
        else
            disp(str);
        end
    end
    function out=disp_error(obj)
        str=sprintf('  Error   %5s = %.3f',obj.Opts.errType,obj.Em);
        if nargout > 0
            out=str;
        else
            disp(str);
        end
    end
    function out=disp_error0(obj)
        str=sprintf('  Error0  %5s = %.3f',obj.Opts.errType,obj.Em0);
        if nargout > 0
            out=str;
        else
            disp(str);
        end
    end
    function init(obj)
        obj.nF=obj.AMA.Index.nF;
        obj.nStim=obj.AMA.Stim.nStim;
        obj.nCtg=obj.AMA.Stim.nCtg;
        obj.ctgInd=obj.AMA.Stim.ctgInd;

        obj.mvncnst=obj.nF*log(2*pi)/2;
        obj.lcorind=num2cell(unique(obj.ctgInd));

        [xcon, ~, ix] = unique(obj.ctgInd);
        nrow = max(ix);
        %ncol = size(,2);
        yrow = (1:obj.nStim)';
        obj.aggIdx = accumarray(ix, yrow, [nrow 1], @(x) {x});
        obj.yagg0 = cell(size(obj.aggIdx));
    end
%- OBJ FUN
    function fun=getFun(obj);
        obj.N.init();
        obj.Opts.init();
        obj.bRan=false;
        switch obj.Opts.ppAlg
        case {'GSS','GMM','FLL','FLLGPU'}
            fun=@(f) obj.objective(f);
        case 'SGD'
            P=obj.getSGDOpts();
            fun=@(f) amaR01objFuncSGD(f,[],P{:});
        case 'GRD'
            P=obj.getGrdOpts();
            fun=@(f) amaR01objFuncFLLgrd(f,[],P{:});
        case 'GRDGPU'
            % TODO
        end
        obj.fun=fun;
    end
    function P=getGrdOpts(obj)
        P={obj.S.s, obj.S.ctgInd, obj.S.X, obj.N.rMax, obj.N.fano, obj.N.var0, obj.Opts.errType};
    end
    function P=getSGDOpts(obj)
        F=obj.FMinOpts;
        P=obj.getSGDOpts;
        P=[obj.N.nF P obj.N.nFset,F.btchSz,F.MaxIter,F.stpSzMax,F.stpSzMin,F.stpSzEta];
    end
%- MAIN
    function Em = objective(obj,f)
        obj.unpack(f);
        if ~obj.bRan
            obj.init();
            obj.bRan=true;
        end
        obj.N.respond();
        obj.posterior();
        [Em,E]=obj.error();
        if isempty(obj.Em0)
            obj.Em0=Em;
        end
    end
    function [Em,E]=test(obj,S)
        %if nargin >= 2 && ~isempty(S)
        %    obj.parse_stim(S);
        %end
        obj.N.respond();
        obj.posterior();
        [Em,E]=obj.error();
    end
    function XHat=estimate(obj,S);
        if nargin > 2 && isempty(S)
            obj.parse_stim(S);
        end
        obj.N.respond();
        obj.posterior();
        XHat=obj.getXHat();

    end
end
methods(Hidden)
%- PACK
    function unpack(obj,fIn)
        [obj.f,obj.w]=obj.Index.unpack(fIn);
    end
%- POSTERIOR
    function posterior(obj)
        switch obj.Opts.ppAlg
        case 'GSS'
            obj.posterior_GSS();
        case 'GMM'
            [obj.PCor,obj.PAll,obj.LCor,obj.LAll] = AMAengineGMM(obj.R,obj.Rm,obj.lambda,obj.S.ctgInd,obj.S.X,obj.NcmpMax); % XXX
        case 'FLL'
                [PCor,PAll]=        AMAengine(obj.N.R,obj.N.Rm,obj.N.nvar, obj.S.ctgInd);
        case 'FLLGPU'
            try
                [obj.PCor,obj.PAll]=AMAengineGPU_comp(obj.N.R,obj.N.Rm,obj.N.nvar, obj.S.ctgInd,max(obj.S.ctgInd));
            catch
                [obj.PCor,obj.PAll]=        AMAengine(obj.N.R,obj.N.Rm,obj.N.nvar,obj.S.ctgInd);
                if ~obj.bErrored
                    obj.disp_gpuError();
                    obj.bErrored=true;
                end
            end
        end
    end
    function dsp_gpuError(obj)
        disp(['amaR01objFuncFLL: WARNING! bGPU = 1, but GPU is not being utilized!!! Quit and type gpuDeviceReset() at command line!']);
    end
    function posterior_GSS(obj)
        % ctgInd [ nStim x 1  ]
        % Rm     [ nStim x nF ]
        % LAll   [ nStim x nCtg ]

        Nr=obj.AMA.Nrn;

        Nr.getLambda();
        L=Nr.Lambda;
        Rm=Nr.Rm;
        R=Nr.R;
        %func= @(x) obj.post_fun(x,R,L);

        % NOTE assuming Rm/stim indeces are presorted

        Xwg=obj.yagg0;
        for iy = 1:length(Xwg)
            Xwg{iy} = obj.post_fun(Rm(obj.aggIdx{iy},:),R,L);
        end
        obj.LAll = exp(cat(2, Xwg{:}) - obj.mvncnst);
        lcor=cellfun(@(i,j) obj.LAll(i,j),obj.aggIdx,obj.lcorind,'UniformOutput',false);
        obj.LCor=cat(1,lcor{:});

        %obj.LAll = cat(2, Xwg{isrt});
        %Xwg = Xwg(isrt,:);
        %size(Xwg)
        Nm=sum(obj.LAll,2);
        obj.PAll = bsxfun(@rdivide,obj.LAll,Nm);
        obj.PCor = bsxfun(@rdivide,obj.LCor,Nm); % @ correct category
    end
    function posterior_GSSold(obj)
        Nr=obj.AMA.Nrn;

        Nr.getLambda();
        L=Nr.Lambda;
        Rm=Nr.Rm;
        R=Nr.R;
        for c = 1:obj.nCtg
            bInd = obj.ctgInd == c; % XXX SLOW 4

            obj.LAll(:,c)=obj.post_fun(Rm(bInd,:),R,L);

            obj.LCor(bInd,1) = obj.LAll(bInd,c); % XXX SLOW 2
        end
    end

    function LLAll= post_fun(~,rm,R,L)
        % sum(log(diag(Rc))) = log determinant

        m  = size(rm,1);
        if m > 1
            s=1;
        else
            s=0;
        end

        mu = sum(rm,1)/m;
        xc = rm - mu;
        Rc = chol( ((xc' * xc) ./ (m-s)) + L);
        %Rc    = obj.cholcov( ((xc' * xc) ./ (m-s)) + L, m);

        X0   = bsxfun(@minus,R, mu);
        LLAll = -0.5 * sum((X0/Rc).^2, 2)  - sum(log(diag(Rc)));
    end
    function T=cholcov(obj,Sigma,m)
        % Test for square, symmetric
        wassparse = issparse(Sigma);
        tol = 10*eps(max(abs(diag(Sigma))));
        if all(all(abs(Sigma - Sigma') < m*tol))
            [T,p] = chol(Sigma);
            if p > 0
                T = zeros(0,'like',Sigma);
            end
        else
            T = zeros(0,'like',Sigma);
        end

        if wassparse
            T = sparse(T);
        end
    end
%- ERROR
    function [Em,E]=error(obj)
        obj.getXHat();
        obj.getXDev();
        obj.getE();
        %obj.regularize();
        if nargout > 0
            Em=obj.Em;
            E=obj.E;
        end
    end
    function Xhat=getXHat(obj)
        switch obj.Opts.estMeth
        case 0
            % MODE
            % TODO
        case 1
            % MEDIAN
            obj.XHat  = zeros(size(obj.PAll,1),1);
            for i = 1:size(obj.PAll,1)
                obj.XHat(i,1) = interp1(cumsum(obj.PAll(i,:)),obj.S.X,0.5,'linear',min(X));
            end
        case 2
            % MEAN
            obj.XHat =  obj.PAll*obj.S.X;
        case 3
            % CIRC
            obj.XHat=circ_mean(obj.S.X,obj.PAll');
        end

        if nargout > 0
            XHat=obj.Xhat;
        end
    end
    function getXDev(obj)
        switch obj.Opts.devMeth
        case 1
            % SUBTRACTIVE
            obj.XDev=obj.XHat-obj.S.XCtg;
        case 2
            % DIVISIVE
            obj.XDev=obj.XHat/obj.S.XCtg;
        case 3
            % CIRC
            obj.XDev=circ_dist(obj.S.XCtg,obj.XHat);
        end
    end
    function getE(obj)
        switch obj.Opts.errMeth
        case -1
            if obj.Opts.expn==0
                % TODO
            else
                E=aexp( obj.XDev, obj.Opts.bAbs+1, obj.Opts.expn);
            end
        case 0
            % L0
            % TODO
        case 1
            % KL-L
            E=-log(obj.LCor);
        case 2
            % KL
            %E=-log(obj.PCor);
            E=-log(obj.PCor);
        case 3
            % MLE
            E=obj.LCor;
        case 4
            % MAP
            E=obj.Pfor;
        end
        obj.E=E .* obj.S.w;
        obj.Em = mean(obj.E);

    end
    function regularize(obj)
        if ~isempty(obj.Opts.regularization);
            %obj.Opts.rlambda=0.1;
            obj.Em=obj.Em + obj.Opts.rlabmda * norm(obj.N.f(:), obj.Opts.regularization);
        end
    end
end
end
