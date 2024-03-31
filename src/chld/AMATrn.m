classdef AMATrn < handle & AMAChld
properties
    Opts

    f00
    f0
    f % XXX MOVE?
    fIn

    c00 % fully raw init
    c0  % init with params
    c
    cIn

    w0
    w % XXX MOVE?
    wIn

    % OUTPUTS
    lambda
    grad
    output
    exitflag
    error
    errorType

    minTimeSec
    fitDate

    UB
    LB

    lb
    ub

    con

% MAKE READ ONLY
end
properties(Access=protected)
    bCmplx

    bAnalytic
    bCmplxParam
    bSplitParam

    nDim
    nPix
    nSplit
    nF


    Prg
    Stim
    Index
    Obj
    ObjOpts

    bInit=false % ???
    lastI
end
properties(Hidden)
    AMADeps={'Nrn','Prg','Index','TrnOpts','Stim','Obj','ObjOpts'}
    lastOpts
end
methods(Static)
end
methods
%- INIT
    function obj=AMATrn(varargin)
        if nargin < 1
            return
        end
    end
    function out=get.Prg(obj)
        out=obj.AMA.Prg;
    end
    function out=get.Index(obj)
        out=obj.AMA.Index;
    end
    function out=get.Stim(obj)
        out=obj.AMA.Stim;
    end
    function out=get.Opts(obj)
        out=obj.AMA.Opts.Trn;
    end
    function out=get.ObjOpts(obj)
        out=obj.AMA.Opts.Obj;
    end
    function out=get.bCmplx(obj)
        out=obj.AMA.Nrn.fCmplx > 0;
    end
    function out=get.bAnalytic(obj)
        out=obj.AMA.Opts.Trn.bAnalytic;
    end
    function out=get.bCmplxParam(obj)
        out=obj.AMA.Opts.Trn.bCmplxParam;
    end
    function out=get.bSplitParam(obj)
        %out=obj.AMA.Opts.Trn.bCmplxParam;
        out=obj.AMA.Opts.Trn.bSplitParam;
    end
    function out=get.nDim(obj)
        out=obj.AMA.Stim.nDim;
    end
    function out=get.nPix(obj)
        out=obj.AMA.Stim.nPix;
    end
    function out=get.nSplit(obj)
        out=obj.AMA.Stim.nSplit;
    end
    function out=get.nF(obj)
        out=obj.AMA.Nrn.nF;
    end
    function out=get.lb(obj)
        out=obj.AMA.Index.lb;
    end
    function out=get.ub(obj)
        out=obj.AMA.Index.ub;
    end
    function assert(obj)
        % CHECK DIMENSIONS
        %assert(isempty(obj.f0) || size(obj.f0,2) == obj.nF,...
        %       'size of f0 dim 2 (%d) does not match number of requested filters (%d).',size(f0,2),nF);
    end
    function parseTrn(obj,varargin)
        obj.lastOpts=obj.Opts.copy([],'AMATrnOpts');
        %if obj.Opts.bParsed
            obj.Opts.parseTrn(varargin{:});
        %else
        %    obj.Opts.parse(varargin{:});
        %    obj.Opts.init();
        %end
        obj.Opts.AMA=obj.AMA;
    end
    function init(obj)
        % TODO add if nF is different
        obj.Opts.init();

        n=obj.Prg.nIter;
        obj.lambda=cell(n,1);
        obj.grad=cell(n,1);
        obj.output=cell(n,1);
        obj.exitflag=ones(n,1)*-3;
        obj.error=zeros(n,1);
        obj.minTimeSec=zeros(n,1);
        obj.fitDate=cell(n,1);


        obj.init_f();
        obj.init_w();
        obj.init_con();
        obj.init_bounds(obj.Opts.fSz);
        obj.bInit=true;
    end
    function [con,lb,ub]=retCons(obj)
        con=obj.con;
        lb=obj.LB;
        ub=obj.UB;
    end
    function checkAndApply(obj)
        obj.Opts.nF > size(obj.fIn,2);
        bApply=~obj.bInit || (~obj.Opts.bContinue && ~obj.Opts.bRecurse && ~obj.Opts.bRecurseW) || obj.Opts.nF > size(obj.fIn,2);
        % RESTORE
        if obj.Opts.bContinue && ~obj.bInit
            obj.Opts=obj.lastOpts;
        end
        % APPLY
        if bApply
            obj.init();
        end
    end
    function setSplitFixInds(obj)
        % fix positive inds for secondary splits
        inds=obj.AMA.Index.getSplitFixInds();
        obj.fIn(inds,:)=0;
        obj.f0(inds,:)=0;

        % unfix parameter ind
        inds=obj.AMA.Index.getSplitParamInds();
        obj.fIn(inds,:)=1;
        obj.f0(inds,:)=1;
    end
    function setAnalyticInds(obj)
        inds=obj.AMA.Index.getAnalyticFixInds();
        obj.fIn(inds,:)=0;
        obj.f0(inds,:)=0;
        if obj.bCmplx
            obj.cIn(inds,:)=0;
            obj.c0(inds,:)=0;
        end
    end
    function setCmplxInds(obj)
        inds=obj.AMA.Index.getCmplxFixInds();
        obj.cIn(inds,:)=0;
        obj.c0(inds,:)=0;

        inds=obj.AMA.Index.getCmplxParamInds();
        obj.cIn(inds,:)=1;
        obj.c0(inds,:)=1;
    end
    function f0=initRandomF(obj)
        f0=AMANrn.initRandomFilters([],obj.Stim.nPix,obj.Opts.nF,obj.Opts.rndSds(1)); % NOTE
    end
    function c0=initRandomC(obj)
        c0=AMANrn.initRandomFilters([],obj.Stim.nPix,obj.Opts.nF,obj.Opts.rndSds(2)); % NOTE
    end
    function loadInF(obj,f)
        if ~isreal(f)
            obj.fIn=real(f);
            obj.cIn=imag(f);
        else
            obj.fIn=f;
        end
        %obj.parameterize();
    end
    function init_f(obj)
        if isempty(obj.fIn)
            nf=0;
        else
            nf=size(obj.fIn,2);
        end
        bDoF=obj.Opts.nF ~= nf;
        if bDoF
            obj.f00 = obj.initRandomF;

            if obj.bCmplx
                obj.c00 = obj.initRandomC();
            end

            % EXTEND
            if obj.Opts.nF > nf
                obj.fIn=[obj.fIn obj.f00(:,nf+1:end)];
                if obj.bCmplx
                    obj.cIn=[obj.cIn obj.c00(:,nf+1:end)];
                end
            % SHRINK
            elseif obj.Opts.nF < nf

                obj.fIn(:,(obj.nF+1):end)=[];
                if obj.bCmplx
                    obj.cIn(:,(obj.nF+1):end)=[];
                end
            end
            obj.f0=obj.f00;
            obj.c0=obj.c00;

            obj.parameterize();

        end

        if obj.Opts.bNewF
            obj.fIn=obj.f0;

            % CMPLX
            if obj.bCmplx
                obj.cIn=obj.c0;
            end
        end
        %if size(obj.fIn,2)==2
        %    figure(88)
        %    imagesc(obj.fIn)
        %    figure(89)
        %    imagesc(obj.cIn)
        %end
    end
    function parameterize(obj)
        % Split
        if obj.bSplitParam && obj.nSplit > 1
            obj.setSplitFixInds();
        end

        % anlytic
        if obj.bAnalytic
            obj.setAnalyticInds();
        end

        % cmplx
        if obj.bCmplx && obj.bCmplxParam
            obj.setCmplxInds();
        end
    end
    function init_w(obj)
        if isempty(obj.w) || any(obj.w(:) < 0)
            nw=0;
            obj.wIn=[];
        else
            nw=size(obj.fIn,2);
        end
        if obj.Opts.nF > nw
            obj.w0 = AMANrn.initRandomWeights([],obj.Opts.nF,obj.Opts.rndSds(3)); % NOTE
            obj.wIn =[obj.wIn obj.w0(nw+1:end,nw+1:end)];
        elseif obj.Opts.nF < nw
            obj.w0 = AMANrn.initRandomWeights([],obj.Opts.nF,obj.Opts.rndSds(3)); % NOTE
        end
        if obj.Opts.bNewW
            obj.wIn=obj.w0;
        end
    end
    function [lb,ub]=init_bounds(obj,fsz)
        lb= ones(fsz)*-1.5;
        ub= ones(fsz)* 1.5;
        if size(fsz,3) > 1
            ub(:,:,2)=1;
            lb(:,:,2)=0;
            inds=obj.nF+1:obj.Stim.nPix;
            ub(inds,:,2)=0;
            lb(inds,:,2)=0;
        end
        obj.LB=lb;
        obj.UB=ub;
    end
    function con=init_con(obj)
        if obj.Opts.bMeanCon
            con=@(rfs) AMATrn.fNMCon(rfs,obj);
        else
            con=@(rfs) AMATrn.fNCon(rfs,obj);
        end
        obj.con=con;
    end
%- PACK
    function pack(obj,I,Prg,varargin)
        % I [f, fval,exitflag,outpt,lambda,grad]

        aInds=obj.Prg.fainds;

        [ft,wt,ct]=obj.Index.separate(varargin{1});
        if obj.AMA.Index.bF
            obj.fIn(:,aInds)=ft;
        end
        if obj.AMA.Index.bW
            obj.wIn(:,aInds)=wt;
        end
        if obj.AMA.Index.bC
            obj.cIn(:,aInds)=ct;
        end
        %obj.fIn(:,obj.Prg.fainds,:)=obj.Index.repack(varargin{1});

        FOut=varargin{1};
        FOut(FOut < obj.AMA.Opts.FMin.TolX)=0;

        [f,w]=obj.Index.unpack(varargin{1});
        if isempty(obj.f)
            obj.f=zeros(size(obj.f00));
        end
        obj.f(:,aInds)=f;
        if ~isempty(w)
            obj.w(:,obj.Prg.wainds)=w;
        end

        % MOVE TO OPTIM ???
        % output may just have everything XXX
        obj.error(I,1) =varargin{2}; %fval

        obj.exitflag(I)=varargin{3};
        obj.output{I}  =varargin{4};
        obj.lambda{I,1}=varargin{5};
        obj.grad{I}    =varargin{6};
        obj.errorType=obj.ObjOpts.errType;

        if ~isempty(Prg)
            obj.minTimeSec(I) = Prg.toc;
            obj.fitDate{I}    = Prg.date;
        else
            obj.minTimeSec(I) = -1;
            obj.fitDate(I)    = -1;
        end
        obj.lastI=I;

        %    % TODO
        %if strcmp(obj.amaType,'SGD')
        %    [obj.Trn.f(:,inds), obj.Ebtch(:,:,i)] = fun(fIndex.F);
        %    return
        %end

    end
end
methods(Static)
    function [c,ceq,cGrd,ceqGrd] = fNMCon(f,obj)
    % with mean==0
        [f,w]=obj.Index.unpack(f);
        f=abs(f);

        c = [];
        % FILTER VECTOR MAGNITUDE (i.e. NORM) MUST EQUAL 1.0
        ceq(1,:) = sum(f.^2)-1;
        ceq(2,:) = mean(f);


        % IF GRADIENT IS PASSED
        if nargout > 2
            nF = size(f,2);
            K = kron(eye(nF),2.*f);
            cGrd = [];
            ceqGrd = K(:,1:(nF+1):end); % GRADIENT OF CONSTRAINT FUNCTION
        end
    end
    function [ciq,ceq,ciqGrd,ceqGrd] = fNCon(f,obj)
        [F,w,C]=obj.Index.unpack(f);
        bCmplx=~isempty(C);
        if bCmplx
            cc=0;
        else
            cc=C;
        end

        % FILTER VECTOR MAGNITUDE (i.e. NORM) MUST EQUAL 1.0
        nF=obj.Index.nf;

        %c=[];
        %ceq(1,:) = sum(f.^2)-(1+bCmplx);

        c   = [];
        ceq = [];

        n=obj.AMA.Stim.PszRC_each;
        nS=obj.AMA.Stim.nSplit;
        N=n*nS;

        I=nS;

        %c1=C(i1,:);
        %c2=C(i2,:);

        %I=(1+bCmplx)*obj.nSplit;
        %I=(1+bCmplx);

        ciqF=zeros(I,nF);
        ciqC=zeros(I,nF);
        for i = 1:I;
            si=obj.AMA.Stim.getSplitInds(i);
            ciqF(i,:)=sum(F(si,:).^2)-1;
            ciqC(i,:)=sum(C(si,:).^2)-1;
        end
        ciq=[ciqF; ciqC];

        % IF GRADIENT IS PASSED
        if nargout > 2
            ceqGrd=[];
            ciqGrd=[];

            g=I.*([F; C]);
            if nF == 1
                ceqGrd=g;
                return
            end
            K = kron(eye(nF), g);
            ciqGrd = K(:,1:(nF+1):end); % GRADIENT OF CONSTRAINT FUNCTION

            %cGrd = [];
            %ceqGrd = K(:,1:(nF+1):end); % GRADIENT OF CONSTRAINT FUNCTION

            %
            %size(cGrd)
        end
    end
end
end
