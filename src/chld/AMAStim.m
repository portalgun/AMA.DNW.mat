classdef AMAStim < handle & AMAChld
properties
    alias
    fname
end
% TODO MAKE READ ONLY
properties
    trainORtest=''
end
properties
    units
    measure
    totT
    totF
    XT
    XF

    % STIM
    S
    As

    % importance weights
    w

    XCtg
    ctgInd
    X

    bBino

    % WHETHER TO
    bMeanSubtract
    bPCATrunc
    bNormalize
    bContrast
    bWhiten
    bWindow
    bFourier

    % WHETHER HAS
    bMeanSubtracted
    bPCATrunced
    bWhitened
    bNormalized
    bContrasted
    bWindowed
    bFouriered

    nPix
    nStim
    PszRC
    PszRC_each
    nDim
    nCtg
    nSplit

    bLoaded=false
    Other
end
properties(Hidden)

    SRaw
    SS
    SF

    bNormalizeS  %whether to fourier
    bFourierS  %whether to fourier

    bPreMeanSubtracted
    bPrePCATrunced
    bPreNormalized
    bPreContrasted
    bPreWhitened
    bPreWindowed
    bPreFouriered

    bParsed=false;

    sTrnAlias
    sTrnFName
    sTstAlias
    sTstFName
end
properties(Hidden)
    allInds
    AMADeps={''}
    hprps={'sTrnAlias','sTrnFName','sTstAlias','sTstFName','bParsed','SS','SF'}
    nprps={'S','alias','fname'}
end
methods(Static)
%- GETP
    function P=getP()
        P=[ ...
            AMAStim.getP_in();
            AMAStim.getP_large();
            AMAStim.getP_base();
            AMAStim.getP_tt();
            AMAStim.getP_do();
            AMAStim.getP_done();
        ];
    end
    function P=getP_in()
        P={...
            'fname','','';
            'alias','','';
        };
    end
    function P=getP_large()
        P={...
            'S',[],'';
            'w',[],'';
            'ctgInd',[],'';
            'X',[],'';
        };
    end
    function P=getP_base()
        P={...
            'bBino',[],'';
            ...
            'measure','','';
            'units','','';
            'totT',1,'';
            ...
            'nDim',[],'';
            'PszRC',[],'';
            'nStim',[],'';
            'nSplit',[],'';
            'nPix',[],'';
            ...
        };
    end
    function P=getP_tt()
        P={ ...
            'sTrnAlias','','';
            'sTstAlias','','';
            'sTrnFName','','';
            'sTstFName','','';
        };
    end
    function P=getP_do()
        P={ ...
            'bMeanSubtract',false,'';
            'bPCATrunc',false,'';
            'bNormalizeS',false,'';
            'bContrast',false,'';
            'bWhiten',false,'';
            'bWindow',false,'';
            'bFourierS',false,'';
        };
    end
    function P=getP_done()
        P={ ...
            'bMeanSubtracted',[],'';
            'bPCATrunced',false,'';
            'bContrasted',false,'';
            'bNormalized',false,'';
            'bWhitened',false,'';
            'bWindowed',[],'';
            'bFouriered',[],'';
        };
    end
%-
    function dire=getDir()
        dire=[Env.var('DATA') 'AMA' filesep 'stim' filesep];
    end
    function out=notTorT(trainORtest)
        if strcmp(trainORtest,'test')
            out='train';
        elseif strcmp(trainORtest,'train')
            out='test';
        end
    end
    function out=getFName(aliasORfname,trainORtest)
        a=[];
        f=[];
        af=aliasORfname;
        tt=trainORtest;
        ntt=AMAStim.notTorT(tt);
        u='_';

        if ~endsWith(af,'.mat')
            mm='.mat';
        else
            mm='';
        end

        % FILE w/ directory
        f0=[af mm];;
        f1=[af u tt mm];
        f2=[af u ntt mm];
        if Fil.is(f1)
            out=f1;
            return
        elseif Fil.is(f0);
            out=f0;
            return
        elseif Fil.is(f2)
            out=f2;
            Error.warnSoft(['Only ' ntt ' stim when requested ' tt]);
            return
        elseif contains(aliasORfname,filesep)
            error('File %s does not exist',aliasORfname);
        end

        % FILE
        dire=AMAStim.getDir();
        f0=[dire af mm];
        f1=[dire af u tt mm];
        f2=[dire af u ntt mm];
        if Fil.is(f1)
            out=f1;
            return
        elseif Fil.is(f0);
            out=f0;
            return
        elseif Fil.is(f2)
            out=f2;
            Error.warnSoft(['Only ' ntt ' stim when requested ' tt]);
            return
        elseif contains(aliasORfname,filesep)
            error('File %s does not exist',aliasORfname);
        end

        % ALIAS
        % LATER do db lookup?

    end
    function obj=Load(aliasORfname,trainORtest)
        obj=AMAStim();
        obj.load_(aliasORfname,trainORtest);
    end
end
methods
    function P=done2do(obj)
        P=obj.getPSave();
        fldstt=AMAStim.getP_done();
        inds=ismember(P(:,1),fldstt(:,1));
        P(inds,1)=regexprep(P(inds,1),'ed$','');
        P{startsWith(P(:,1),'bNormaliz'),1}='bNormalize';
    end
    function P=getPSave(obj)
        P=[ ...
            AMAStim.getP_base();
            AMAStim.getP_tt();
            AMAStim.getP_done();
        ];
    end
    function obj=AMAStim(varargin)
        if nargin < 1
            return
        end
    end
    function init(obj)
        if ~obj.bLoaded
            obj.load();
        end
    end
    function obj=load_old(obj,aliasORfname)
        % ARGS IN
        if nargin < 2
            AF=[];
        else
            AF=aliasORfname;
        end
        if nargin < 3 || isempty(trainORtest)
            if ~isempty(obj.trainORtest)
                trainORtest=obj.trainORtest;
            else
                trainORtest=[];
            end
        end

        % FNAME & ALIAS
        if isempty(AF)
            if ~isempty(obj.fname)
                fname=obj.fname;
            elseif  ~isempty(obj.alias)
                alias=obj.alias;
            else
                error('No filename provided')
            end
        else
            fname=AMAStim.getFName(AF);
            alias=AMAStim.fname2alias(obj.fname);
        end
        if isempty(fname)
            error('No stimulus filename or alias provided')
        end
        % SET
        obj.fname=fname;
        obj.alias=alias;
    end
    function load(obj)
        if isempty(obj.fname) && ~isempty(obj.alias)
            fname=AMAStim.getFName(obj.alias,obj.trainORtest);

            switch obj.trainORtest
            case 'test'
                obj.sTstFName=fname;
            case 'train'
                obj.sTrnFName=fname;
            end

        elseif isempty(obj.fname) && isempty(obj.alias)
            error('No fname or alias provided')
        end

        S=load(obj.fname);
        clear cl;
        flds=fieldnames(S);
        % Unnest
        while numel(flds)==1 && isstruct(S.(flds{1}))
            S=S.(flds{1});
            flds=fieldnames(S);
        end
        if isfield(S,'s')
            stm=S.s;
            S=rmfield(S,'s');
        elseif isfield(S,'S')
            stm=S.S;
            S=rmfield(S,'s');
        end
        obj.parse(S,stm);
    end
    function parse(obj,S,stm)

        % RM FIELDS
        rmflds={'trainORtest','fname'};
        flds=fieldnames(S);
        rmflds(~ismember(rmflds,flds))=[];
        S=rmfield(S,rmflds);

        % OPTS
        [~,obj.Other]=Args.simpleLooseExclusive(obj,AMAStim.getP(),S);


        % STIM
        % XXX verify stim are sorted by ind
        obj.ctgInd=Vec.col(obj.ctgInd);
        obj.nStim=numel(obj.ctgInd);
        if size(stm,1)==obj.nStim && size(stm,2)~=obj.nStim
            stm=stm';
        end
        obj.nPix=size(stm,1);

        if isempty(obj.bFouriered)
            obj.bFouriered=~all(isreal(stm));
        end

        % ASSIGN STIM
        obj.SRaw=stm;
        if obj.bFouriered
            obj.SF=stm;
        else
            obj.SS=stm;
        end

        if isempty(obj.bMeanSubtracted) && ~obj.bFouriered
            obj.bMeanSubtracted=all(mean(obj.SS,1) < 10^-8);
        end
        if isempty(obj.bWindowed) && ~obj.bFouriered
            obj.bWindowed=all([obj.SS(1,:) obj.SS(end,:)] < 10^-6);
        end

        % PRE PROCESSED STATUS
        flds=obj.getProcFlds(true);
        for f = 1:length(flds)
            fld=['b' flds{f}];
            pfld=['bPre' flds{f}];
            obj.(pfld)=obj.(fld);
        end

        % DEFAULTS
        if isempty(obj.X) && ~isempty(obj.ctgInd)
            obj.X=unique(obj.ctgInd);
        end

        if isempty(obj.w)
            obj.w=1;
        end
        if isempty(obj.PszRC)
            obj.PszRC=size(obj.S,1);
        end
        if isempty(obj.nSplit)
            if ~isempty(obj.bBino);
                obj.nSplit=obj.bBino+1;
            elseif isfield(S,'measure')
                switch S.measure
                case {'disparity','Disparity'}
                    obj.nSplit=2;
                otherwise
                    obj.nSplit=1;
                end
            else
                obj.nSplit=1;
            end
        end


        obj.nCtg=numel(obj.X);
        obj.nDim=sum(obj.PszRC > 1);
        if ~isempty(obj.X)
            obj.XCtg=obj.X(obj.ctgInd);
        end
        if obj.nDim==1
            obj.PszRC_each=[obj.PszRC(1)./obj.nSplit 1];
        elseif obj.nDim==2
            XXX
        end
        obj.allInds=obj.getSplitInds(1:obj.nSplit);



        % XXX get dimensions smartly
        obj.getXTF();

        obj.proc();

        obj.bParsed=true;
        obj.bLoaded=true;
    end
    function getXTF(obj)
        if obj.bFouriered && isempty(obj.totF)
            obj.totF=1;
        elseif ~obj.bFouriered && isempty(obj.totT)
            obj.totT=1;
        end

        NX=obj.PszRC_each(1:obj.nDim);
        if obj.bFourierS
            obj.totT=NX./obj.totF;
        else
            obj.totF=NX./obj.totT;
        end
        for i = 1:obj.nDim
            [obj.XT(i,:)]=Ft.getXT(NX(i),obj.totT(i),[],[],true,true);

            [obj.XF(i,:)]=Ft.getXT(NX(i),obj.totF(i),[],[],true,true);
        end
    end
%- PROCESS
    function reset(obj)
        if obj.bFouriered
            obj.SF=obj.SRaw;
        else
            obj.SS=obj.SRaw;
        end

        flds=obj.getProcFlds(true);
        for f = 1:length(flds)
            fld=['b' flds{f}];
            pfld=['bPre' flds{f}];
            obj.(fld)=obj.(pfld);
        end

    end
    function out=bTransform(obj)
        flds=obj.getProcFlds(false);
        eflds=obj.getProcFlds(true);
        out=false;
        for i = 1:length(flds)
            %fld=flds{i};
            fld=['b' flds{i}];
            efld=eflds{i};
            if strcmp(fld,'bFourier')
                continue
            end

            if obj.(fld) && ~obj.(efld)
                out=true;
                return
            end

        end
    end
    function proc(obj)
        obj.reset();

        if obj.bFouriered && obj.bTransform()
            obj.ifft_();
        end
        if obj.bMeanSubtract && ~obj.bMeanSubtracted
            obj.meanSubtract();
        end
        if obj.bPCATrunc && ~obj.bPCATrunced
            obj.PCATrunc();
        end
        if obj.bContrast && ~obj.bContrasted
            obj.contrast();
        end
        if obj.bNormalize && ~obj.bNormalized
            obj.normalize();
        end
        if obj.bWhiten && ~obj.bWhitened
            obj.whiten();
        end
        if obj.bWindow && ~obj.bWindowed
            obj.window();
        end
        if obj.bFourierS && ~obj.bFouriered
            obj.fft_();
        end
    end
%-PROCESS ROUTINES
    function N=contrast(obj)
        obj.SS=bsxfun(@rdivide,obj.SS,mean(obj.SS));
        obj.bContrasted=true;
        obj.bContrast=false;
    end
    function normalize(obj)
        obj.SS=bsxfun(@rdivide,obj.SS,sqrt(sum(obj.SS.^2)));
        obj.bNormalized=true;
        obj.bNormalizeS=false;
    end
    function meanSubtract(obj)
        obj.SS=bsxfun(@minus,obj.SS,mean(obj.SS));
        obj.bMeanSubtracted=true;
        obj.bMeanSubtract=false;
    end
    function PCATrunc(obj)
        % XXX
        return
        lambda = pca(obj.SS');
        obj.bPCATrunced=true;
        obj.bPCATrunc=false;
    end
    function whiten(obj)
        X=obj.SS';
        %mu = mean(X);
        %X = bsxfun(@minus, X, mu);
        A = X'*X;
        [V,D,~] = svd(A);
        W = sqrt(size(X,1)-1)*V*sqrtm(inv(D + eye(size(D))*epsilon))*V';
        Xwh = X*W;
        %invMat = pinv(whMat);

        obj.SS=XWh';
        obj.bWhitened=true;
        obj.bWhiten=false;
    end
    function window(obj)
        if obj.bWindowed
            return
        end
        env=(Envl.getEnv('cos',obj.XT,[],0.5)-.5)*2;
        env=Vec.col(env);
        if obj.nSplit==2
            env=[env; env];
        end

        %XXX handle pre fourier

        ss=obj.SS;
        ss=bsxfun(@times,ss,env);
        obj.SS=ss;

        obj.bWindowed=true;
        obj.bWindow=false;
    end
    function fft_(obj,SS)
        % S    [ nPix  x nStm ]
        [obj.SF,obj.As]=obj.fft(obj.SS);
        obj.bFouriered=true;
        obj.bFourier=false;
    end
    function [SF,As]=fft(obj,SS)
        if obj.nDim == 1
            [SF,As]=obj.fft1D(SS);
        else
            [SF,As]=obj.fft2D(SS);
        end
    end
    function SF=ifft(obj,SS)
        if obj.nDim == 1
            SF=obj.ifft1D(SS);
        else
            SF=obj.ifft2D(SS);
        end
    end
    function [SF,As]=fft1D(obj,SS)
        % SS [nPix nDim]
        N=obj.PszRC_each(1);
        SF=zeros(size(SS));
        As=zeros(size(SS));

        for i = 1:obj.nSplit
            inds=(1:N)+(N*(i-1));
            SF(inds,:)=fftshift( fft( ifftshift( SS(inds,:),1) ,[],1) ,1);  % HERE FFT
            As(inds,:)=abs(SF(inds,:));
        end
    end
    function SS=ifft1D(obj,SF)
        % SS [nPix nDim]
        N=obj.PszRC_each(1);
        SS=zeros(size(SF));

        for i = 1:obj.nSplit
            inds=(1:N)+(N*(i-1));
            SS(inds,:)=ifftshift( ifft( ifftshift( SF(inds,:),1) ,[],1) ,1);  % HERE FFT
        end
    end
    function SF=fft2D(obj,SF)
        % XXX
    end
    function SF=ifft2D(obj,SF)
        % XXX
    end
%- ASSERT
    function assert(obj)
        assert(prod(obj.PszRC)==obj.nPix,...
            'fSize product =%d must be equal to nPix=%d',prod(obj.PszRC),obj.nPix);
        assert(obj.nFset >= 1 && obj.nFset <= obj.nF,...
               'nFset=%d. Must be greater than 0 and less than nF= ',obj.nF);
        assert(obj.nF < obj.nPix,...
               'amaR01: WARNING! more filters requested nF=%d than stimulus dimensionality size(s,1)=%d', obj.nF, obj.nPix);
    end
%- SET/GET
    function out=get.bNormalize(obj)
        out=obj.bNormalizeS;
    end
    function out=get.bFourier(obj)
        out=obj.bFourierS;
    end
    function set.bFourier(obj,val)
        obj.bFourierS=val;
    end
    function set.bNormalize(obj,val)
        obj.bNormalizeS=val;
    end
    function out=get.alias(obj)
        switch obj.trainORtest
        case 'test'
            out=obj.sTstAlias;
        case 'train'
            out=obj.sTrnAlias;
        otherwise
            out='';
        end
    end
    %function out=get.fname(obj)
    %    out=obj.sfname;
    %end
    %function set.fname(obj,val)
    %    obj.sfname=val;
    %end
    function out=get.fname(obj)
        switch obj.trainORtest
        case 'test'
            out=obj.sTstFName;
        case 'train'
            out=obj.sTrnFName;
        otherwise
            out=obj.sTrnFName;
        end
    end
    function set.fname(obj,val)
        switch obj.trainORtest
        case 'test'
            obj.sTstFName=val;
        case 'train'
            obj.sTrnFName=val;
        otherwise
            obj.sTrnFName=val;
        end
    end
    function set.alias(obj,val)
        if isempty(val)
            return
        end
        switch obj.trainORtest
        case 'test'
            out=obj.sTstAlias;
        case 'train'
            out=obj.sTrnAlias;
        otherwise
            error(['invalid alias ' val]);
        end
        obj.bLoaded=false;
    end
    function set.bFourierS(obj,in)
        obj.bFourierS=in;
        if ~obj.bParsed
            return
        end
        if in && isempty(obj.SF)
            obj.fft_();
        end
    end
    function out=get.S(obj)
        if obj.bFourierS
            out=obj.SF;
        else
            out=obj.SS;
        end
    end
%- INDS %
    function split(obj,splitI)
    end
    function inds=getSplitInds(obj,splitI,n)
        N=obj.PszRC_each(1);
        if nargin < 2 || isempty(splitI)
            splitI=1;
        end
        if nargin < 3 || isempty(n)
            n=N;
        end
        ii=(1:n);
        inds=zeros(n,length(splitI));
        for i =1:length(splitI)
            inds(:,i)=ii+(N*(splitI(i)-1));
        end
    end
    function out=getDCInds(obj,splitI)
        if nargin < 2 || isempty(splitI)
            splitI=1:obj.nSplit;
        end
        N=obj.PszRC_each(1);

        out=floor(N/2+1)+(N.*(splitI-1));
    end
    function inds=getNegInds(obj,splitI,bZero,bPad)
        if nargin < 2 || isempty(splitI)
            splitI=1:obj.nSplit;
        end
        if nargin < 3
            bZero=false;
        end
        if nargin < 3
            bPad=false;
        end
        N=obj.PszRC_each(1);
        if bPad || mod(N,2)~=0
            ss=1;
        else
            ss=2;
        end

        dcI=obj.getDCInds(splitI)-~bZero;
        bI=ss+N.*(splitI-1);
        inds=zeros(floor(N/2+bZero),length(dcI));
        for i = 1:length(dcI);
            inds(:,i)=bI(i):dcI(i);
        end
    end
    function inds=getPosInds(obj,splitI,bZero)
        if nargin < 2 || isempty(splitI)
            splitI=1:obj.nSplit;
        end
        if nargin < 3
            bZero=true;
        end
        N=obj.PszRC_each(1);

        dcI=obj.getDCInds(splitI)-bZero;
        bI=N.*(splitI);
        inds=zeros(bI(1)-dcI(1)+1,length(dcI));
        for i = 1:length(dcI);
            inds(:,i)=dcI(i):bI(i);
        end
    end
%- OPERATIONS
    function reshape(obj,sz)
        % XXX
    end
    function out=dot(obj,A)
        out=dot(obj.S,A);
    end
    function out=mtimes(obj,A,bSplit)
        %s       [nPix  x nStm]
        %f       [nPix  x nF]
        %r       [nStim x nF]
        % TODO correct dims from start

        if nargin < 3 || isempty(bSplit)
            bSplit=false;
        end
        if bSplit
            out=zeros(obj.nStim,size(A,2)*obj.nSplit);
            % TODO MAKE NOT LOOP
            for i = 1:obj.nSplit
                inds=obj.allInds(:,i);
                out(:,i)=obj.S(inds,:)'*A(inds,:);
            end
        else
            out=obj.S'*A;
        end
    end
%% UTIL
    function flds=getProcFlds(obj,bEd)
        if nargin < 2
            bEd=false;
        end
        flds={'MeanSubtract','Normalize','Whiten','Window','Fourier'};
        if ~bEd
            return
        end
        for i = 1:length(flds)
            fld=flds{i};
            if endsWith(fld,'e')
                flds{i}=[fld 'd'];
            else
                flds{i}=[fld 'ed'];
            end
        end
    end
    function lbl=getxlabel(obj)
        if isempty(obj.measure)
            X='X';
        else
            X=obj.measure;
        end
        if isempty(obj.units)
            units='';
        else
            units=[' (' obj.units ')'];

        end
        lbl=[X units];
    end
    function varargout=index(obj,ind,bSplit,bFourier)
        if nargin < 2 || isempty(ind)
            ind=1;
        end
        if nargin < 3 || isempty(bSplit)
            bSplit=obj.nSplit > 1;
        end
        if nargin < 4 || isempty(bFourier)
            bFourier=obj.bFourierS;
        end
        if bFourier
            s=obj.SF(:,ind);
        else
            s=obj.SS(:,ind);
        end
        if obj.nDim==1
            s=s';
        end
        if bSplit
            n=obj.PszRC_each;
            varargout{1}=s(:,1:n);
            varargout{2}=s(:,(n+1):end);
        else
            varargout{1}=s;
        end
    end
end
end
