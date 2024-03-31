classdef AMAIndex < handle & AMAChld
properties
    nF
    nPix
    sz
    nf

    f0
    c0
    w0

    fFix
    fLrn

    cFix
    cLrn

    wFix
    wLrn

    skipFFix
    skipWFix

    % W OPTS
    bSym
    bCmpl
    bEye

    bF
    bC
    bW

    % PACKED
    F
    %V

    FInd
    CInd
    WInd

    rmInds
end
% XXX Make read only
properties(Hidden)
    UB
    LB

    lb
    ub

    bAnalytic
    bCmplxParam
    bSplitParam

    bCmplx
    nSplit
    PszRC_each
end
properties(Access=protected)
    TrnOpts
    Trn
    Prg
    Nrn
end
properties(Hidden)
    AMADeps={'TrnOpts','Trn','Prg','Nrn','Stim'}
end
methods
    function obj=AMAIndex(varargin)
        if nargin < 1
            return
        end
        obj.pack(varargin{:});
    end
    function out=get.Trn(obj)
        out=obj.AMA.Trn;
    end
    function out=get.TrnOpts(obj)
        out=obj.AMA.Opts.Trn;
    end
    function out=get.Prg(obj)
        out=obj.AMA.Prg;
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
        out=obj.AMA.Opts.Trn.bCmplxParam;
    end

    function out=get.nSplit(obj)
        out=obj.AMA.Stim.nSplit;
    end
    function out=get.PszRC_each(obj)
        out=obj.AMA.Stim.PszRC_each;
    end

    function out=get.LB(obj)
        out=obj.AMA.Trn.LB;
    end
    function out=get.UB(obj)
        out=obj.AMA.Trn.UB;
    end

    function pack(obj,Trn,Prg)
        f=Trn.fIn;
        c=Trn.cIn;
        %figure(808)
        %imagesc(Trn.f0)
        %dk

        w=Trn.wIn;
        fIndsLrn=Prg.finds;
        fIndsFix=Prg.ffinds;
        wIndsLrn=Prg.winds;
        wIndsFix=Prg.wfinds;

        obj.nPix=size(f,1);
        obj.nF=size(f,2);

        fInds=sort([fIndsLrn fIndsFix]);
        wInds=sort([wIndsLrn wIndsFix]);
        obj.bW=~isempty(wInds);
        obj.bF=~isempty(fInds);
        obj.bC=~isempty(fInds) && ~isempty(c);
        if ~obj.bW & ~obj.bF
            error('Nothing set')
        elseif ~obj.bW || ~obj.bF
            ;
        elseif ~isequal(fInds,wInds)
            error('All indeces need accounting for');
        end

        I=0;
        if obj.bF
            I=I+1;
            obj.FInd=I;
        end
        if obj.bC
            I=I+1;
            obj.CInd=I;
        end
        if obj.bW
            I=I+1;
            obj.WInd=I;
        end

        % rminds are filters indeces not yet being used
        obj.rmInds=~ismember(1:obj.nF,fInds);

        obj.F=[];
        obj.init_f(f,fIndsLrn,fIndsFix);
        obj.init_w(w,wIndsLrn,wIndsFix);
        obj.init_c(c,fIndsFix);


        obj.sz=size(obj.F);
        obj.nf=size(obj.F,2);

    end
    function [f0,lb,ub]=retF0(obj)
        obj.modifyBounds();
        lb=obj.lb(:);
        ub=obj.ub(:);
        f0=obj.F(:);
    end
    function [lb,ub]=modifyBounds(obj)
        i1=1:size(obj.F,1);
        i2=1:size(obj.F,2);
        i3=1:size(obj.F,3);
        lb=obj.LB(:,i2,i3);
        ub=obj.UB(:,i2,i3);

        for k = i2

            bFix=~isempty(obj.AMA.Prg.ffinds) && ismember(k,obj.AMA.Prg.ffinds);

            c=1;
            if bFix
                I=i1; % FIX FILTER
            else
                I=obj.fFix(:,k); % FIX SOME INDECES
            end
            %m=Vec.col(find(I));
            %l=repmat(k,size(m,1),1);
            %ind=sub2ind(size(obj.f0),m,l);

            lb(I,k,c)=obj.f0(I,k);
            ub(I,k,c)=obj.f0(I,k);

            if obj.bC
                c=c+1;
                if bFix
                    I=i1;
                else
                    I=obj.cFix(:,k);
                end

                lb(I,k,c)=obj.c0(I,k);
                ub(I,k,c)=obj.c0(I,k);
                if ~bFix
                    lb(~I,k,c)=-1;
                    ub(~I,k,c)= 1;
                    lb(I,k,c)= 0;
                    ub(I,k,c)= 0;
                end
            end
            if obj.bW
                bFix=~isempty(obj.AMA.Prg.wfinds) && ismember(k,obj.AMA.Prg.wfinds);
                c=c+1;
                if obj.AMA.Prg.wfinds==k
                    I=i1;
                else
                    I=obj.wFix(:,k);
                end

                lb(I,k,c)=obj.w0(I);
                ub(I,k,c)=obj.w0(I);
            end
        end
        obj.lb=lb;
        obj.ub=ub;

        obj.checkBounds();
    end
    function checkBounds(obj)
        %tol=10^-7;
        tol=obj.AMA.Opts.FMin.TolCon;
        q=((obj.lb - obj.F) > tol) | ((obj.ub-obj.F) < -tol);

        if any(q)
            obj.plotBounds();
            obj.plotBad();

            disp([obj.lb(q)  obj.ub(q) obj.F(q)])
            error('f0 not within bounds')
        end

    end
    function plotBounds(obj)
        nP=size(obj.lb,3);

        Fig.new('AMA Bounds and F0');
        for p = 1:nP
            obj.plotBoundsFun_(obj.lb,p,1);
            obj.plotBoundsFun_(obj.ub,p,2);
            obj.plotBoundsFun_(obj.F,p,3);
        end
    end
    function plotBad(obj)
        nP=size(obj.lb,3);
        tol=obj.AMA.Opts.FMin.TolCon;
        bLB=((obj.lb-obj.F) >  tol);
        bUB=((obj.ub-obj.F) < -tol);

        Z=zeros(size(obj.F));
        lb=Z;
        ub=Z;
        lb(bLB)=1;
        ub(bUB)=1;

        Fig.new('AMA Bad Bounds and F0');
        for p = 1:nP
            obj.plotBoundsFun_(lb,p,1,2);
            obj.plotBoundsFun_(ub,p,2,2);
            %obj.plotBoundsFun_(obj.F,p,3);
        end
    end
    function plotBoundsFun_(obj,b,p,n, nn)
        if nargin < 5 || isempty(nn)
            nn=3;
        end
        nF=size(obj.lb,2);
        nP=size(obj.lb,3);
        switch p
            case obj.FInd
                titl='f';
            case obj.CInd;
                titl='c';
            case obj.WInd;
                titl='c';
            otherwise
                titl='';
        end
        switch n
        case 1
            lbl='lb';
        case 2
            lbl='ub';
        case 3
            lbl='F0';
        end

        subPlot([nn nP],n,p);
        hold off;
        imagesc(indexFun(b,p));
        caxis([-1 1]);
        axis square;
        if n==1
            title(titl);
        end
        if p == 1
            ylabel(lbl);
        end
        function out=indexFun(in,p)
            out=[];
            nF=size(in,2);
            for j = 1:nF
            for i = 1:obj.nSplit
                inds=obj.AMA.Stim.getSplitInds(i);
                out=[out in(inds,j,p)];
            end
            end
            %out=out';
        end
    end
    function inds=getAnalyticFixInds(obj)
        inds=obj.AMA.Stim.getNegInds([],true,true);
    end
    function inds=getSplitMainInds(obj)
        inds=obj.AMA.Stim.getPosInds(1,false);
    end
    function inds=getSplitFixInds(obj)
        inds=obj.AMA.Stim.getPosInds(2:obj.nSplit,false);
    end
    function inds=getSplitParamInds(obj)
        inds=obj.AMA.Stim.getDCInds(2:obj.nSplit)+1;
    end
    function inds=getCmplxFixInds(obj)
        inds=obj.AMA.Stim.getPosInds([],true);
    end
    function inds=getCmplxParamInds(obj)
        inds=obj.AMA.Stim.getDCInds(1)+(1:obj.nSplit);
    end
    function fixSplit(obj)
        % fix positive inds for secondary splits
        inds=obj.getSplitFixInds();
        obj.fFix(inds,:)=true;

        % unfix parameter ind
        inds=obj.getSplitParamInds();
        obj.fFix(inds,:)=false;
    end
    function fixAnalytic(obj)
        inds=obj.getAnalyticFixInds();
        obj.fFix(inds,:)=true;
    end
    function fixCmplx(obj)
        % fix
        inds=obj.getCmplxFixInds();
        obj.cFix(inds,:)=true;

        inds=obj.getCmplxParamInds();
        obj.cFix(inds,:)=false;
    end
    function constrainW(obj)
        % SYMMETRICIZE W
        obj.bSym=false; %TODO
        obj.bCmpl=false; % TODO
        if obj.bSym | obj.bCmpl
            obj.wLrn=obj.wLrn & triu(true(size(obj.wLrn)));
            obj.wFix=obj.wFix | tril(true(size(obj.wLrn)));
        end

        % DIAGONALIZE W
        obj.bEye=false; % TODO
        if obj.bEye
            obj.wLrn=Arr.onEye(obj.wLrn,false);
            obj.wFix=Arr.onEye(obj.wLrn,true);
        end
    end
    function init_f(obj,f,fIndsLrn,fIndsFix)
        zf=false(obj.nPix,obj.nF,1);

        % F Lrn
        obj.fLrn=zf;
        obj.fLrn(:,fIndsLrn)=true;

        % F Fix
        obj.fFix=zf;
        obj.fFix(:,fIndsFix)=true;
        obj.fFix(:,obj.rmInds)=[];


        % Secondary Splits
        if obj.bSplitParam && obj.nSplit > 1
            obj.fixSplit();
        end

        % Analytic (fix negative Inds to 0)
        if obj.bAnalytic
            obj.fixAnalytic();
        end
        %obj.fFix(:,fIndsFix)=true;
        %obj.fFix(obj.rmInds,:)=[];

        obj.skipFFix=~any(any(obj.fFix));

        f(:,obj.rmInds)=[];
        obj.f0=f;
        obj.F=cat(3,obj.F,f);
    end
    function init_c(obj,c,fIndsFix)
        if ~obj.bC
            return
        end

        obj.cFix=obj.fFix;
        obj.cLrn=obj.cLrn;

        % Parameterize Complex
        if obj.bC && obj.bCmplxParam
            obj.fixCmplx();
        end

        obj.cFix(:,fIndsFix)=true;

        c(:,obj.rmInds)=[];
        obj.c0=c;
        obj.F=cat(3,obj.F,c);
    end
    function init_w(obj,w,wIndsLrn,wIndsFix)
        zw=false(obj.nF,obj.nF,1);

        % W LRN
        obj.wLrn=zw;
        obj.wLrn(wIndsLrn,:)=true;


        % W FIX
        obj.wFix=zw;
        obj.wFix(wIndsFix,:)=true;
        obj.wFix(:,obj.rmInds)=[];
        obj.wFix(obj.rmInds,:)=[];
        obj.skipWFix=~any(any(obj.wFix));

        if ~obj.bW
            return
        end

        % CONSTRAIN W
        obj.constrainW();

        w(:,obj.rmInds)=[];
        w(obj.rmInds,:)=[];

        nFNew=numel(fInds);
        nz=obj.nPix-nFNew;
        w=[w; zeros(nz,nFNew)];
        obj.w0=w(1:size(w,2),:);

        obj.F=cat(3,obj.F,w);
    end
    function [fIn]=repack(obj,fOut)
        [f,w,c]=obj.separate(fOut);
        fIn=f;
        for i = 2:3
            if     obj.bC && obj.CInd==i
                fIn=cat(3,fIn,c);
            elseif obj.bW && obj.WInd==i
                fIn=cat(3,fIn,w);
            end
        end
    end
    function [f,w,c]=separate(obj,fIn)
        fIn=reshape(fIn,obj.sz);

        %- UNPACK PARTS
        % f
        f=fIn(:,:,obj.FInd);

        % w
        if obj.bW
            w=fIn(:,:,obj.WInd);
        else
            w=[];
        end

        % c
        if obj.bC
            c=fIn(:,:,obj.CInd);
        else
            c=[];
        end
    end
    function [f,w,c]=unpack(obj,fIn)
        c=[];
        [fraw,w,c]=obj.separate(fIn);

        % SET FIXED f
        if ~obj.skipFFix && obj.bF
            fraw(obj.fFix)=obj.f0(obj.fFix); % this is correct
            if obj.bC
                c(obj.cFix)=obj.c0(obj.cFix);
            end
        end

        if obj.bF & obj.bSplitParam && obj.nSplit > 1
            % SPLIT PARAMATERIZATION
            fixI  =obj.getSplitFixInds();
            paramI=obj.getSplitParamInds();
            mainI =obj.getSplitMainInds();
            ff=fraw;
            fraw(fixI,:) = bsxfun(@times, fraw(mainI,:), fraw(paramI,:));

            %subplot(1,2,1);
            %imagesc(ff);

            %subplot(1,2,2);
            %imagesc(f);
            %drawnow

        end

        if obj.bF & obj.bC && obj.bAnalytic
            % COMPLEX PARAMETERIZATION
            if obj.bCmplxParam
                fc=fraw;
                m=c(~obj.cFix);
                N=obj.PszRC_each(1);
                for i = 1:obj.nSplit
                    sinds=(1:N)+(N*(i-1));
                    fc(sinds,:)=fc(sinds,:)*m(i);
                end
                if nargout > 2
                    f=fraw;
                    c=fc;
                else
                    f=complex(fraw,fc);
                end
            else
                if nargout > 2
                    f=fraw;
                else
                    f=complex(fraw,c);
                end
            end
        else
            f=fraw;
        end



        if ~obj.bW; return; end

        %- UNPACK WEIGHTS
        if obj.bW
            w=fraw(1:obj.nF,:,obj.WInd);
        else
            w=1;
        end

        % SET FIXED w
        if ~obj.skipWFix && obj.bW
            w(obj.wFix)=obj.w0(obj.wFix);
        end

        % APPLY W SYMMETRY CONSTRAINTS
        if obj.bCmpl
            w=Arr.triuOnTrilCmpl(w);
        elseif obj.bSym
            w=Arr.triuOnTril(w);
        end
        if obj.bEye
            w=Arr.onEye(w,0);
        end

    end
    function [f,w]=unpackLrn2(obj)

        nFLrn=sum(obj.fLrn);
        nWLrn=sum(obj.wLrn);
        f=obj.f0;
        w=obj.w0;
        if nFLrn > 0
            f(obj.fLrn)=obj.F(1:fLrn);
        end
        if nWLrn > 0
            w(obj.wLrn)=obj.F(fLrn+1:end);
        end
        if obj.bCmpl
            w=Arr.triuOnTrilCmpl(w);
        elseif t0.bSym
            w=Arr.triuOnTril(w);
        end
        if nargout < 1
            obj.f=f;
            obj.w=w;
        end
    end
end
end
