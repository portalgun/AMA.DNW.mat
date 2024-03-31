classdef AMA_Stats < handle
properties
end
methods(Static)
end
methods
    function wav(obj)
        i=1;
        [Lx,Rx]=obj.get_Pctrs();
        while i <= obj.nStim
            [LCor,R]=obj.splitS(i);
            C=wcoherence(LCor,R);

            figure(203)

            subPlot([1 2],1,1);
            hold off;
            plot(LCor,'r'); hold on
            plot(R,'k');
            plot([Rx(i) Rx(i)],[-0.5 0.5],'r');
            plot([Lx(i) Lx(i)],[-0.5 0.5],'k');
            axis square;
            ylim([-0.5 .5]);
            Axis.format('X','Lum');

            subPlot([1 2],1,2);
            hold off;
            imagesc(C);
            hold on;

            y=size(C,1)+1;
            plot([Rx(i) Rx(i)],[0 y],'r');
            plot([Lx(i) Lx(i)],[0 y],'k');
            Axis.format('x','f');
            colormap bone;
            axis square;
            colorbar;
            caxis([0 1]);
            i=plotFlipper(i);
        end
    end
    function P=powerSpec(obj)
        [Lx,Rx]=obj.get_Pctrs();
        nPix=obj.S.nPix;
        wdw=cosWdw.new([nPix/2 nPix/2],[nPix/2 nPix/2],.5,[1 1]);
        ind=ceil(nPix/2);
        wdw=wdw(:,ind);
        for i = 1:obj.nStim
            [l,r]=obj.splitS(i);

            %C(:,i)=mscohere(LCor,R,wdw);
            LCor(i,:)=pwelch(l);
            R(i,:)=pwelch(r);
            X(i,:)=cpsd(l,r);
        end
        P.Co=real(X);
        P.Q=imag(X);

        P.LCor=LCor;
        P.R=R;
        P.X=X;

    end
    function P2=powerSpec2(obj)
        P=obj.A.P;
        flds=fieldnames(P);
        for i = 1:length(flds)
            fld=flds{i};
            P2.(['m' fld])=mean(P.(fld),2);
            P2.(['v' fld])=var(P.(fld),[],2);
        end

        P2.MSC=abs(P.X).^2./(P.LCor.*P.R);
        P2.mMSC=(mean(P.Co,2).^2+mean(P.Q,2).^2)./mean(P.LCor.*P.R);
        P2.mX=mean(P.Co,2)+mean(P.Q,2);
        P2.sMSC=sum(P2.MSC,2);
        P2.vMSC=var(P2.MSC,[],2);

        %P2.MSCL=P2.MSC.*(P.LCor+P.R);
    end
    function dsp2=getDisparityContrast(obj,inds)
        nPix=obj.S.nPix;
        wdw=cosWdw.new([nPix/2 nPix/2],[nPix/2 nPix/2],.5,[1 1]);
        wdw=wdw(nPix/4,:);
        dsp2=zeros(obj.S.nStim,1);
        if nargin < 2
            inds=1:obj.S.nStim;
        end
        for I = 1:numel(inds)
            i=inds(I);
            dsp2(i)=XYZ.disparity_contrast(obj.S.Ixyz(:,:,:,i),obj.S.PszRC,'dim',1,'bConv',false,'winPosXYZm',[0 0 1],'W',wdw);
        end
        obj.A.dsp2=dsp2;
    end
    function bi=getBiDiff(obj)
        bi=zeros(obj.S.nStim,1);
        nPix=obj.S.nPix;
        wdw=cosWdw.new([nPix/2 nPix/2],[nPix/2 nPix/2],.5,[1 1]);
        wdw=wdw(nPix/4,:);
        wv=wdw(:);
        ws=sum(wv);
        for i = 1:obj.S.nStim
            [LCor,R]=obj.splitS(i);
            bi(i)=sqrt(mean( (R(:)-LCor(:)).^2./wv)/ws);
        end
        obj.A.bi=bi;
    end
    function [OUT]=XCorr(obj)
        [Lx,Rx]=obj.get_Pctrs();

        m=max(obj.S.PszRC)-1;
        n=length(Lx);

        [LCor,R]=obj.splitS(1);

        AL=zeros(m,n);
        AR=zeros(m,n);
        C=zeros(m,n);
        Crnk=zeros(m,n);
        [~,lag]=xcorr(LCor,R);
        for i = 1:n
            [LCor,R]=obj.splitS(i);

            C(:,i)=xcorr(LCor,R);

            %[~,Crnk(:,i)]=sort((C(:,i)),'descend');
            [~,Crnk(:,i)]=sort(abs(C(:,i)),'descend');

            %[~,Crnk(:,i)]=ismember(C(:,i),sort(C(:,i),'descend'));
            lagRnk(:,i)=lag(Crnk(:,i));
            %AL(:,i)=xcorr(LCor,LCor);
            %AR(:,i)=xcorr(R,R);
            %Ch(:,i)=mscohere(LCor,R);
            %Xh(:,i)=cpsd(LCor,R);
            %PL(:,i)=pwelch(LCor);
            %PR(:,i)=pwelch(R);
        end

        for i = 2:n
            [LCor,R]=obj.splitS(i);
        end

        [Lx,Rx]=obj.get_Pctrs();
        ELag=round(Rx-Lx);

        % Max
        mxs=max(C,[],1);
        inds=find_fun(C,mxs);
        MaxOLag=lag(inds)';
        MaxDiffLag=(MaxOLag-ELag);


        % Mean
        mn=mean(C,1);
        v=var(C,[],1)';

        MeanOCorr=mean(C,1)';

        %- CORR
        OUT.mxInds=inds;
        OUT.C=C;
        W=1;
        mx=max(abs(lag));
        Em=[0; repelem((1:mx)',2,1)];
        for i = 1:obj.nStim
            lg(:,i)=abs(abs(lagRnk(:,i))-(Em+ELag(i))).*W;
        end
        lagErr=mean(lg,1);

        OUT.lag=lag;
        OUT.CRnk=Crnk;
        OUT.lagRnk=lagRnk;
        OUT.lagErr=lagErr;


        OUT.Lx=Lx;
        OUT.Rx=Rx;
        OUT.ELag=ELag;
        OUT.MeanOCorr=MeanOCorr;
        OUT.MaxOLag=MaxOLag;
        OUT.MaxDiffLag=MaxDiffLag;

        function [inds,bMult]=find_fun(C,in)
            n=length(in);
            inds=zeros(n,1);
            bMult=false(n,1);
            for i = 1:n
                ind=find(C(:,i)==in(i));
                if numel(ind)==1
                    inds(i)=ind;
                else
                    bMult(i)=true;
                    inds(i)=ind(1);
                end
            end
        end
    end
    function Lu=LumStats(obj)
        Lu=struct();
        n=obj.nStim;
        Lu.vl=zeros(n,1);
        Lu.vr=zeros(n,1);
        Lu.vb=zeros(n,1);
        Lu.vB=zeros(n,1);

        Lu.mr=zeros(n,1);
        Lu.ml=zeros(n,1);
        Lu.mb=zeros(n,1);
        for i = 1:n
            [LCor,R]=obj.splitS(i);
            Lu.vl(i)=var(LCor);
            Lu.vr(i)=var(R);
            Lu.vb(i)=var([LCor; R]);
            Lu.vB(i)=mean([var(LCor) var(R)],2);
            Lu.mB(i)=mean([mean(LCor) mean(R)],2);


            Lu.ml(i)=mean(LCor);
            Lu.mr(i)=mean(R);
            Lu.mb(i)=mean([LCor; R]);
        end
        Lu.diffS=Lu.vl-Lu.vr;
        Lu.diffM=Lu.ml-Lu.mr;
    end
    function clearStat(obj,names)
        if ~iscell(names)
            names={names};
        end
        for i = 1:length(names)
            obj.A.(names{i})=[];
        end
    end
    function getStats(obj,bForce)
        if nargin < 2
            bForce=false;
        end
        obj.clearStat('P2');
        flds={'X','XCorr';
              'LCor','LumStats';
              'P','powerSpec';
              'P2','powerSpec2';
        };
        for i = 1:size(flds,1)
            fld=flds{i,1};
            meth=flds{i,2};
            if ~isfield(obj.A,fld) || isempty(obj.A.(fld)) || bForce
                obj.A.(fld)=obj.(meth);
            end
        end
        [~,obj.A.Em.MAP]=obj.getError('MAP');
        [~,obj.A.Em.MES,obj.A.Em.XHat]=obj.getError('MES');
        [~,~,~,obj.A.Em.bC]=obj.getNCorrect();
        obj.A.Em.bGd=obj.A.Em.bC | abs(obj.A.Em.MES) <= 1 | obj.A.Em.MAP <= 1;
        obj.A.Em.X=obj.S.X(obj.S.ctgInd);
    end
    function errAnalysis(obj)
        obj.getStats();
        A=obj.A;
        Em=obj.A.Em;
        X=A.X;
        P=A.P2;
        LCor=A.LCor;


        obj.fig('X Err New');
        plot3(Em.X,A.bi,Em.MAP,'.');

        %%- X ERR
        obj.fig('X Err');
        hold off;

        %- LINE
        m=diff(obj.S.X(2:-1:1));
        m=-1;
        xint=mean(Em.X);
        xint=-8.1;
        b=-m*xint;
        l=@(x) m*x + b;
        lx=[min(Em.X),max(Em.X)];
        y=l(lx);

        %- CURVE
        for i = 1:length(obj.S.X)
            X=obj.S.X(i);
            inds=obj.S.ctgInd==i;
            [c(i,:),e(i,:)]=hist(Em.MES(inds),25);
            mx(i)=max(c(i,:));
            my(i)=e(i,find(c(i,:)==mx(i),1,'first'));
        end

        ct=obj.S.ctgInd;
        XC=obj.S.X(ct);
        inds1=((Em.MES > l(XC) & XC < xint) | ...
               (Em.MES < l(XC) & XC > xint) ...
        );
        inds2=((Em.MES < my(ct)' & XC < xint) | ...
               (Em.MES > my(ct)' & XC > xint) ...
        );
        obj.A.inds1=inds1;
        obj.A.inds2=inds2;

        plot(lx,y,'r','LineWidth',3); hold on
        plot(obj.S.X,my,'r','LineWidth',3); hold on
        obj.plotGd(Em.X,Em.MES,[],inds1|inds2);
        plot(lx,[0 0],'k','LineWidth',3);
        Axis.format('X','MES Erro');


        %obj.fig('X Err New');
        %hold off;
        %obj.plotGd(Em.X,P.m,(Em.MES));
        %Axis.format3('X','C','MES Erro');


        %%- X ERR Hist
        obj.fig('X Err Hist');
        hold off;
        for i = 1:length(obj.S.X)
            x=obj.S.X(i)*ones(25,1);
            plot3(x,e(i,:),c(i,:),'LineWidth',2); hold on
        end
        plot3(obj.S.X,my,mx,'k','LineWidth',2);
        plot3(lx,y,zeros(size(y)),'k','LineWidth',3); hold on
        plot3([min(obj.S.X),max(obj.S.X)],[0 0],[0 0],'k','LineWidth',3);
        Axis.format3('X','Err','Count');
    end
    function plotGd(obj,f1,f2,f3,bC);
        if nargin < 5 || isempty(bC)
            bC=obj.A.Em.bC;
        end
        if nargin < 4 || isempty(f3)
            plot(f1(~bC), f2(~bC),'k.'); hold on;
            plot(f1( bC), f2( bC),'b.'); hold on;
        elseif nargin >=4
            hold off;
            if sum(~bC) > 0
                plot3(f1(~bC), f2(~bC), f3(~bC),'k.'); hold on;
            end
            if sum(bC) > 0
                plot3(f1( bC), f2( bC), f3( bC),'b.'); hold on;
            end
        end
    end
    function LumAn(obj,Lu)
        [~,errMAP]=obj.getError('MAP');
        [~,errMES]=obj.getError('MES');
        errMAP=errMAP';

        %plot(errMES,abs(Lu.vl-Lu.vr),'.');
        %plot(errMES,abs(Lu.ml-Lu.mr),'.');
        %plot3(errMAP,abs(Lu.ml-Lu.mr),abs(Lu.vl-Lu.vr),'.');
        plot(errMAP,abs(Lu.vl-Lu.vr),'.');
        plot(errMAP,abs(Lu.ml-Lu.mr),'.');
        %plot(errMAP,Lu.vl,'.');
        %plot(errMES,Lu.vl,'.');
        %plot(errMAP,Lu.vb,'.');
        plot(errMAP,Lu.vB,'.');
        plot(errMAP,abs(Lu.mB),'.');
        plot3((Lu.mB),Lu.vB,errMAP,'.');
        %hist(abs(Lu.mB)./errMAP);
        Axis.format('','');
        zlabel('error');
        %Axis.format('Est. Error','Max Lag Error');
    end
    function XCorrAn(obj,C,lag,Lu)
        if nargin < 2
            [C,lag]=obj.XCorr;
        end
        if nargin < 3
            Lu=obj.LumStats();
        end
        [~,errMAP]=obj.getError('MAP');
        [~,errMES]=obj.getError('MES');

        [Lx,Rx]=obj.get_Pctrs();


        N=3;
        c=0;

        figure(989)
        c=c+1;
        subPlot([1 N],1,c);
        plot(errMES,(MaxDiffLag),'.');
        Axis.format('Est. Error','Max Lag Error');

        c=c+1;
        subPlot([1 N],1,c);
        plot(errMAP,MaxDiffLag,'.');
        Axis.format('Est. Error','Mean Corr');

        c=c+1;
        subPlot([1 N],1,c);
        plot(errMAP,abs(MeanOCorr),'.');
        Axis.format('MAP Est. Error','Mean Corr');
        xlim([1 6]);

        figure(990)
        %plot3(errMAP,abs(MeanOCorr),abs(Lu.vl-Lu.vr),'.');
        %Axis.format('MAP Est. Error','Mean Corr','lumDiff');
        %plot(abs(MeanOCorr),abs(Lu.ml-Lu.mr),'.');


        %plot(VarOLag,err,'.');
        %Axis.format('Var Lag','Est. Error');
        %Axis.format('OLag Var','Est. Error');
    end
    function plotWCohere(obj,f)
        if nargin < 2 || isempty(f)
            fld=obj.getOFld;
            f=obj.(fld).f;
        end
        [c,cp,X,Y]=RNorm.getCohereW(f);

        obj.fig('Cohere W');
        imagesc(c);
        axis image;
        Axis.format('','outputs','inputs');
        colorbar;
        colormap hot;

        %obj.fig('Cohere Sign');
        %imagesc(c>c');
        %axis image;
        %Axis.format('','outputs','inputs');
        %colorbar;
        %colormap hot;

        %obj.fig('Cohere mag');
        %imagesc(abs(c));
        %axis image;
        %Axis.format('','outputs','inputs');
        %colorbar;
        %colormap hot;

        obj.fig('Cohere WM');
        imagesc(cp);
        axis image;
        Axis.format('','outputs','inputs');
        colorbar;
        colormap hot;


    end
end
end
