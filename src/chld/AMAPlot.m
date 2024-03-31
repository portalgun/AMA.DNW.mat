classdef AMAPlot < handle & AMAChld
properties
    bParentPlot=true
    TSNE=struct();
end
properties(Hidden)
    AMADeps={'Nrn','Stim'}
end
properties(Access=protected)
    Nrn
    Stim
    Prg
end
methods(Static)
end
methods
    function out=get.Nrn(obj)
        out=obj.AMA.Nrn;
    end
    function out=get.Prg(obj)
        out=obj.AMA.Prg;
    end
    function out=get.Stim(obj)
        out=obj.AMA.Stim;
    end
%- RESPONSES
    function R(obj,fPairs,ctg2plt,bCombine)
        if nargin < 2
            fPairs=[];
        end
        if nargin < 3
            ctg2plt=[];
        end
        if nargin < 4
            bCombine=[];
        end
        P=obj.initRParms();
        Fig.newIfBase('AMA R');
        obj.R_(fPairs,ctg2plt,bCombine,[],[],P);
    end
    function RStep(obj,fPairs,ctg2plt,bCombine)
        if nargin < 2
            fPairs=[];
        end
        if nargin < 3
            ctg2plt=[];
        end
        if nargin < 4
            bCombine=[];
        end
        P=obj.initRParms();
        P.pauseT=0.2;
        Fig.newIfBase('AMA R');
        obj.R_(fPairs,ctg2plt,bCombine,[],[],P);
    end
    function tsne(obj,fR,fF,bReplot)
        if nargin < 2
            fR='Rm';
        end
        if nargin < 3
            fF=[];
        end
        if nargin < 4 || isempty(bReplot)
            bReplot=false;
        end
        Fig.newIfBase('AMA tsne');
        obj.tsne_(fF,fR,bReplot);
    end
    function retsne(obj,fR,fF)
        if nargin < 2
            fR='Rm';
        end
        if nargin < 3
            fF=[];
        end
        obj.tsne_(fF,fR,true);
    end
%- FILTERS
    function F(obj,nFset,ctg2plt)
        if nargin < 2
            nFset=[];
        end
        if nargin < 3
            bFFT=[];
        end
        obj.F_(nFset,bFFT,[],[]);
    end
    function FS(obj,nFset,ctg2plt)
        if nargin < 2
            nFset=[];
        end
        Fig.newIfBase('AMA FS');
        obj.F_(nFset,false,[],[]);
    end
    function FF(obj,nFset,ctg2plt)
        if nargin < 2
            nFset=[];
        end
        Fig.newIfBase('AMA FF');
        obj.F_(nFset,true,[],[]);
    end
    function F0(obj,nFset,ctg2plt)
        if nargin < 2
            nFset=[];
        end
        if nargin < 3
            bFFT=[];
        end
        Fig.newIfBase('AMA F0');
        obj.F_(nFset,bFFT,[],0);
    end
    function F00(obj,nFset,ctg2plt)
        if nargin < 2
            nFset=[];
        end
        if nargin < 3
            bFFT=[];
        end
        Fig.newIfBase('AMA F0');
        obj.F_(nFset,bFFT,[],-3);
    end
    function FIn(obj,nFset,ctg2plt)
        if nargin < 2
            nFset=[];
        end
        if nargin < 3
            bFFT=[];
        end
        Fig.newIfBase('AMA FIn');
        obj.F_(nFset,bFFT,[],-4);
    end
    function FTrn(obj)
        Fig.new('AMA-train filters');
        obj.bParentPlot=false;
        obj.F_([],[],[],2);

        %ctg2plt=[];
        %fl=obj.Prg.finds;
        %ff=obj.Prg.ffinds;
        %if numel(f)==1
        %    if isempty(ff)
        %        fPairs=[fl fl+1];
        %    else
        %        fPairs=[fl max(ff)];
        %    end
        %elseif numel(f)==2
        %end


        %Fig.new('AMA-train response');
        %obj.R_(fPairs,[],[],[],[],[],false);
    end
%- STIM
    function SS(obj,ind)
        if nargin < 2 || isempty(ind)
            ind=[];
        end
        obj.stim(ind,false);
    end
    function SF(obj,ind)
        if nargin < 2 || isempty(ind)
            ind=[];
        end
        obj.stim(ind,true);
    end
end
methods(Hidden)
    function [R,Nrn]=getR(obj,fF,fR,bCombine)
        if nargin < 2
            fF=[];
        end
        if nargin < 3
            fR=[];
        end
        if nargin < 4
            bCombine=[];
        end
        [R,Nrn]=obj.AMA.getR(fF,fR,bCombine);
    end
    function f=getF(obj,fF)
        if nargin < 2
            fF=[];
        end
        f=obj.AMA.getF(fF);
    end
    function P=initRParms(obj)
        P=struct('bPlotEllipse',true,'bPltRsp',true,'bPlotMarg',true,'pauseT',0,'axisLims',[]);
    end
    function F_(obj,nFset,bFFT,bCmplx,fF,f)
        % -    - real
        % :    - imaginary
        % red  - L
        % blue - R

        if nargin < 3
            bFFT=[];
        end
        if nargin < 5
            fF=[];
        end



        %plotFilters_Disparity(f,fSmp,[])
        %return
        % F
        if nargin < 6 || isempty(f)
            f=obj.getF(fF);
            bCustom=false;
        else
            bCustom=true;
        end

        if nargin < 2 || isempty(nFset)
            nFset = size(f,2);
        end

        S=obj.Stim;
        bSameSplit=S.nSplit > 1; % XXX

        nPix=S.nPix;
        fSize=S.PszRC;
        nFSplit=S.nSplit;
        b2D=S.nDim > 1;

        bFourier=obj.AMA.Nrn.bFourierF;

        if isempty(bFFT)
            bFourierPlot=bFourier;
        else
            bFourierPlot=bFFT;
        end

        if bCustom
            bIFFT=false;
            bFFT=false;
        elseif bFourier
            xt=S.XF;
            if isempty(bFFT) || bFFT
                bFFT=false;
                bIFFT=false;
            elseif ~bFFT
                bFFT=false;
                bIFFT=true;
            end
        else
            xt=S.XT;
            if isempty(bFFT) || bFFT
                bFFT=true;
                bIFFT=false;
            elseif ~bFFT
                bIFFT=false;
                bIFFT=false;
            end
        end

        if bFourierPlot
            xlbl='Frequency';
            xt=S.XF;
        else
            xlbl='Position';
            xt=S.XT;
        end
        len=max(xt)-min(xt);

        if nargin < 4 || isempty(bCmplx)
            bCmplx=obj.Nrn.fCmplx > 0;
        end

        % FFT
        if bFFT
            f=S.fft(f);
        elseif bIFFT
            f=S.ifft(f);
        end


        % CMPLX
        if ~bCmplx
            %f=abs(f);
        end

        %hold off; plot(real(f(:,1))); hold on;  plot(imag(f(:,1)));
        %dk

        %15, 32

        % ylim
        maxi=max(abs(f(:)));
        rn=2*maxi*.1;
        yl=[-(maxi+rn) maxi+rn];

        mm=6/sqrt(nPix/nFSplit);
        %yl=[-mm mm];

        % xlim
        xl0=max(abs(xt)).*[-1 1];
        xl=xl0;
        nx=numel(xt);
        %if bSameSplit
        %    %xl(2)=xl(2)+len*(nFSplit-1)-1;
        %end

        % SUB-DIMS
        if b2D
            % XXX
            M=nF;
            N=nFset;
        else
            M=ceil(sqrt(nFset));
            N=max([ceil(nFset/M),1]);
            colors=['r','b','k','m'];
        end

        % PLOT
        szs=S.PszRC_each;
        d=szs(1);
        for ii = 1:nFset

            % EYE-SPECIFIC FILTER INDEXES
            % LEFT AND RIGHT EYE FILTERS
            for j = 1:nFSplit
                xti=xt;
                if j == 1
                    inds=1:d;
                else
                    inds=(d+1):size(S.S,1);
                    %if bSameSplit
                    %    xti=xt+len*(j-1);
                    %end
                end
                %inds = (1:d)+(j-1)*d;

                if b2D
                    sp=subPlot([N,M],ii,j);
                else
                    if j==1;
                        sp=subPlot([M,N],ii,1);
                        spPos=sp.Position;
                        delete(sp);
                    end
                end
                if bSameSplit
                    w=spPos(3)/nFSplit;
                    h=spPos(4);
                    x=spPos(1)+w*(j-1);
                    y=spPos(2);
                    axes('Position',[x y w h]);
                end
                if b2D
                    imagesc(reshape(f(inds,ii),fSize));
                    % XXX
                else
                    clr=colors(j);

                    %plot(1:numel(inds),f(inds,ii),[colors(j) '-']);
                    fi=f(inds,ii);
                    if bCmplx
                        plot(xti,real(fi),[clr  '-']); hold on
                        plot(xti,imag(fi),[clr  ':'], 'LineWidth',2);
                        %plot(xt, abs(fi), [clr  '-'],'LineWidth',2);
                    else
                        plot(xti,real(fi),[clr '-']);
                        hold on;
                    end
                    if bFourierPlot
                        plot([0 0],yl,'k:');
                    end
                end
                if bSameSplit
                    format_fun(j==1 & ii==nFset,j==1 && ii==nFset,j==1,j==1);
                end
            end
            if ~bSameSplit
                format_fun();
            end
            %if bSameSplit
            %    xticks(xtiks);
            %    xticklabels(xtikl);
            %end
        end
        function format_fun(bXlbl,bYlbl,bTitle,bYTicks)
            if nargin < 1
                bXlbl=true;
            end
            if nargin < 2
                bYlbl=true;
            end
            if nargin < 3
                bTitle=true;
            end
            if nargin < 4
                bYTicks=true;
            end
            hold off;

            axis square;
            if b2D
                caxis([mini maxi]);
            else
                axis xy;
                if ~bFFT
                    xlim(xl);
                    ylim(yl);
                end
            end
            titl=[];
            xL=[];
            yL=[];
            if bTitle
                titl=['RF=' num2str(ii)];
            end
            if bXlbl
                xL=xlbl;
            end
            if bYlbl
                yL='Weight';
            end
            Axis.format(xL,yL,titl);
            xtls=xticklabels();
            if ~bYTicks
                yticklabels([]);
            end
            if bSameSplit
                obj.xtick_split_fun_(j);
            end
        end
    end
    function R_(obj,fPairs,ctg2plt,bCombine,fF,fR,P, f)
        % fPairs:       filter pairs to plot                 [ nPairs  x    2    ]
        % axisLims:     manual control over axis limits (e.g. rMax*[-1 1 -1 1])

        % COMPUTE FILTER RESPONSES TO TRAINING AND TEST STIMS


        if nargin < 5 || isempty(bCombine)
            bCombine=false;
        end
        if nargin < 5 || isempty(fF)
            fF=[];
        end
        if nargin < 6
            fR=[];
            fR='Rc';
        end


        % NRN
        [RR,Nrn]=obj.getR(fF,fR,bCombine);

        bSplit=Nrn.bSplitR;
        bSplit=true; % XXX note

        X=Vec.row(obj.Stim.X);
        ctgInd=obj.Stim.ctgInd;

        rMax=obj.Nrn.rMax;

        if nargin < 2 || isempty(fPairs)
            if Nrn.fCmplx==3 && size(RR,2)==2
                fPairs=[];
            elseif size(RR,2)>1
                fPairs=[1 2];
            else
                fPairs=[1 1];
            end
        end
        if Nrn.fCmplx==3
            if isempty(fPairs)
                fPairs=1:2;
            else
                A=fPairs;
                B=fPairs;
                A(A > 1)=A(A > 1)+1;
                B(B > 1)=B(B > 1)+2;
                B(1)=B(1)+1;
                fPairs=[A; B]';
                if bSplit;
                    fPairs=[fPairs; fPairs'];
                end
            end
        end
        %fPairs=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4]

        if nargin < 3 || isempty(ctg2plt)
            ctg2plt=unique(ctgInd);
        end


        % XXX extraOpts
        bLegend = false;
        CI      = 90;
        modelCR = 'gaussian';


        % MLE FIT CONDITIONAL RESPONSE DISTIRBUTIONS
        if P.bPlotEllipse
            [MU,COV,SIGMA,KAPPA] = fitCondRespDistribution(modelCR,[],[],RR,ctgInd,rMax);
        end

        % COLRS
        colorss = getColorOrder([],length(ctg2plt));
        colors(ctg2plt(:),:) = colorss;

        %fPairs
        %size(RR)
        bNeg=any(RR(:,unique(fPairs)) < 0);


        f=gcf;
        nP=ceil(size(fPairs,1));
        nC=length(ctg2plt);
        %if nStim.nSplit > 1
        %    getSplitInds()
        %end
        for t = 1:nP
            R =RR(:,fPairs(t,:));
            R1=R(:,1);
            R2=R(:,2);

            set(groot,'CurrentFigure',f.Number);
            hold off;
            subplot(1,nP,t);
            hold off;

            % AXIS LIMS
            if isempty(P.axisLims)
                if t == 1 & bNeg
                    lims=[-1 1 -1 1];
                else
                    lims=[-.1 1 -.1 1];
                end
                P.axisLims = 1.2.*max(max(abs(R))).*lims;
            end

            % PLOT CONDITIONAL DISTRIBUTIONS
            for c = 1:nC;
                ctg=ctg2plt(c);
                clrs=colors(ctg,:);

                ind     = find(ctgInd==ctg); % & sqrt(sum(rTst(:,1:2).^2,2)) < 0.2);
                indRnd  = randsample(1:length(ind),min([1000 numel(ind)]));

                r1=R1(ind(indRnd));
                r2=R2(ind(indRnd));

                % ELLIPSE
                if P.bPlotEllipse
                    he(c) = plotEllipse(MU(ctg,:),COV(:,:,ctg),CI,fPairs(t,:),2,colors(ctg,:));
                    hold on;
                end
                % SCATTER
                if P.bPltRsp
                    ht(c) = plot(r1,r2,'wo','linewidth',.125,'markerface',clrs,'markersize',6);
                    hold on;
                end

                % MARGINALS
                if P.bPlotMarg
                    hm(c)=obj.hist_fun(r1,r2,P.axisLims,0.1,clrs,1,31);
                    hold on;
                end

                % COMBINED MARGINALS
                if P.bPlotMarg && c == nC
                    obj.hist_fun(R1,R2,P.axisLims,0.2,[0 0 0],2.0,100);
                    hold on;
                end

                % FORMAT
                xlbl=['F' num2str(fPairs(t,1)) ' response'];
                ylbl=['F' num2str(fPairs(t,2)) ' response'];
                titl=['X=' num2str(X(ctg),'%.3f')];
                Axis.format(xlbl,ylbl,titl);
                %
                axis(P.axisLims);
                axis square;

                % PAUSE
                if P.pauseT > 0
                    pause(P.pauseT);
                end

            end
            hold off;
        end
        if P.bPlotEllipse
            h=he;
        elseif bPLotRsp
            h=ht;
        end
        hold off;

        % LEGEND
        if bLegend
            set(groot,'CurrentFigure',f.Number);
            lbl=legendLabel('X=',X(ctg2plt),1,4);
            legend(h,lbl,'Location','NorthEast');
        end
        hold off;

    end
    function RMarg(obj,fF,fR,bCombine)
        if nargin < 2
            fF=[];
        end
        if nargin < 3
            fR='Rm';
        end
        if nargin < 4
            bCombine=true;
        end
        % XXX
        %bCombine=true

        Stim=obj.Stim;
        X=Stim.X;
        ctgInd=Stim.ctgInd;
        n=Stim.nSplit;

        if bCombine
            nn=0;
        else
            nn=n-1;
        end

        [R,Nrn]=obj.getR(fF,fR,bCombine);
        Fig.newIfBase('AMA RMarg');
        hold off;

        if bCombine
            N=size(R,2);
        else
            N=size(R,2)/n;
        end

        bInc=false;

        [~,edges]=hist(R(:));
        cnts=zeros(numel(edges),N,nn+1,length(X));

        for o = 1:N
        for q = 0:nn
            s=o+q;
            for x =1:length(X)
                r=R(ctgInd==x,s);
                [cnts(:,o,q+1,x),~]=hist(r,edges);
            end
        end
        end
        yl=[0 max(cnts(:))];

        bSplit=n > 1 && ~bCombine;

        for o = 1:N
            if bSplit
                spPos=obj.sp_pos(N,o);
            end
            for q = 0:nn
                s=o+q;
                if bSplit
                    ax=obj.sp_split_fun(spPos,q+1);
                else
                    subPlot([N,nn+1],o,q+1); hold off;
                end

                for x =1:length(X)
                    r=R(ctgInd==x,s);

                    if bCombine
                        c=cnts(:,o,:,x);
                        c=c(:);
                    else
                        c=cnts(:,o,q+1,x);
                    end

                    plot(edges,c); hold on
                end
                axis square;
                if bSplit
                    if o==1 && q==0
                        titl='L';
                    elseif o==1 && q==1
                        titl='R';
                    else
                        titl='';
                    end
                else
                    titl='';
                end

                Axis.format([],['RF= ' num2str(o)],titl);
                ylim(yl);
                if ~bInc
                    hold on;
                end
            end

            if bInc
                drawnow;
                pause(0.5);
            end
        end
    end
%- TSNE
    function fit_tsne_(obj,fF,fR)
        % t-Distributed Stochastic Neighbor Embedding
        if nargin < 2
            fF=[];
        end
        if nargin < 3 || isempty(fR)
            fR='Rm';
        end
        distM='cosine';
        distM='correlation';
        distM='hamming';
        distM='euclidean';

        R=obj.getR(fF,fR);
        obj.TSNE.(fR)=tsne(R,'Distance',distM,'Standardize',false,'Perplexity',50);
    end
    function tsne_(obj,fF,fR,bReplot)
        if nargin < 2
            fF=[];
        end
        if nargin < 3 || isempty(fR)
            fR='Rm';
        end
        if nargin < 4
            bReplot=false;
        end

        % FIT
        if bReplot || ~isfield(obj.TSNE,fR) || isempty(obj.TSNE.(fR))
            obj.fit_tsne_(fF,fR);
        end

        f=Fig.new(['AMA tsne ' fR]);
        hold off;

        Y=obj.TSNE.(fR);
        hold off;
        gscatter(Y(:,1),Y(:,2),obj.Stim.XCtg,[],[],8);
        Axis.format('','');
        axis square;
        %Fig.suptitle(['t-SNE Responses ' newline obj.get_title()]);
    end
%- STIM
    function stim(obj,ind,bFourier)
        if nargin < 2 || isempty(ind)
            ind=1;
        end
        if nargin < 3
            bFourier=[];
        end

        S=obj.AMA.Stim;
        [s1,s2]=S.index(ind,[],bFourier);

        if bFourier
            x=S.XF;
            s=S.SF;
        else
            x=S.XT;
            s=S.SS;
        end

        c=S.ctgInd(ind);
        xC=S.X(c);
        titl=sprintf('Stim=%d Ctg=%d',ind,c);
        if S.nSplit==2
            if S.nDim==1


                spPos=obj.sp_pos(1,1);
                ax=obj.sp_split_fun(spPos,1);

                plot(x,real(s1),'k'); hold on
                plot(x,imag(s1),'k:'); hold on
                Axis.format();
                %xlabel(S.getxlabel);
                title(titl);
                axis square;
                obj.xtick_split_fun_(j);

                ax=obj.sp_split_fun(spPos,2);

                plot(x,real(s2),'k'); hold on
                plot(x,imag(s1),'k:'); hold on
                Axis.format();
                yticklabels({});
                axis square;
                obj.xtick_split_fun_(j);

                %xlabel(S.getxlabel);
            end
        end
    end
    function spPos=sp_pos(obj,n,i)
        sp=subplot(n,1,i);
        spPos=sp.Position;
        delete(sp);
    end
    function ax=sp_split_fun(obj,spPos,j)
        w=spPos(3)/obj.Stim.nSplit;
        h=spPos(4);
        xp=spPos(1)+w*(j-1);
        y=spPos(2);
        pos=[xp y w h];

        ax=axes('Position',[xp y w h]);
    end
    function xtick_split_fun_(obj,j)
        if any(ismember(xlim,xticks))
            xtl=xticklabels();
            if j==1
                n=length(xtl{end});
                spc=repmat(' ',1,n*2);
                xtl{end}=[xtl{end} spc];
            else
                n=length(xtl{1});
                spc=repmat(' ',1,n*2);
                xtl{1}=[spc xtl{1}];
            end
            xticklabels(xtl);
        end
    end
end
methods(Static)
    function h=hist_fun(r1,r2,axisLims,m,colr,lw,n)
        ll=linspace(axisLims(1),axisLims(2),n);
        [H1,B1]=hist(r1,ll);
        [H2,B2]=hist(r2,ll);
        ly=m.*diff(axisLims(1:2)).*H1./max(H1) + axisLims(1);
          plot(B1,ly,'color',colr,'lineWidth',lw);
        h=plot(ly,B2,'color',colr,'lineWidth',lw);
    end
end
end
