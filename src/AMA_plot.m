classdef AMA_plot < handle
properties
end
methods(Static)
end
methods
%- UTIL
    function fnew=splitF(obj)
        % XXX MOVE
        nPixE=obj.nPix/obj.nFSplit;
        d=floor(nPixE);
        fnew=zeros(d,obj.nFSplit,obj.nF);
        for j = 1:obj.nFSplit
            inds = (1:d)+(j-1)*d;
            fnew(:,j,:)=obj.(fld).f(inds,:,:);
        end
        fnew=fnew(:,:);
    end
    function [LCor,R]=splitS(obj,i)
        s=obj.S.s(:,i);
        e=obj.nPix/obj.nFSplit;
        LCor=s(1:e);
        R=s(e+1:end);
    end
    function f=fig(obj,titl,bNew,WH,bFocus)
        if nargin < 3 || isempty(bNew)
            bNew=true;
        end
        if nargin < 4
            WH=[];
        end
        if nargin < 5
            bFocus=bNew;
        end

        if isempty(obj.bParentPlot)
            obj.bParentPlot=0;
        end
        if ~obj.bParentPlot
            a=strrep(obj.alias,'_',' ');
            T=[titl ' ' a];
            f=Fig.new(T,'bClear',bNew,'bFocus',bFocus);
            if bNew && ~isempty(WH)
                set(f,'Position',[f.Position(1:2) WH]);
            end
        else
            f=gcf;
        end
    end
    function [Lx,Rx]=get_Pctrs(obj)
        pOpts=obj.S.ptchOpts;
        assignin('base','pOpts',pOpts);
        vOpts=View3D.optsFromPtchOpts(pOpts);
        View3D.clearCache();

        X=obj.S.X/60;
        if obj.S.PszRC(2)==1;
            PszRC=[obj.S.PszRC(1) obj.S.PszRC(1)]/2;
        else
            PszRC=obj.S.PszRC;
        end
        Lx=zeros(length(X),1);
        Rx=zeros(length(X),1);
        for i = 1:length(X)
            vOpts.Dsp=X(i);
            V=View3D(vOpts);
            PctrCPs=V.get_patch_CPs(PszRC,[]);
            Lx(i)=PctrCPs{1}(2);
            Rx(i)=PctrCPs{2}(2);
        end
        Lx=Lx(obj.S.ctgInd);
        Rx=Rx(obj.S.ctgInd);
    end
    function out=get_title(obj)
        if isempty(obj.alias)
            alias=obj.genAlias();
        else
            alias=obj.alias;
        end
        if strcmp(obj.trainORtest,'test')
            s=obj.sAalias;
        else
            s=obj.stAlias;
        end
        out=[strrep(alias,'_','-') newline   strrep(s,'_','-') ' ' upper(obj.trainORtest)];
    end
%- MAIN
    function testExp(obj)
        %N=linspace(.8,3,10);
        N=1:.2:4;
        bStart=true;
        for i = N
            warning('off');
            normType=['c__atten_R' num2str(i)];
            obj.plotREachCtg(normType,0,bStart);
            bStart=false;
            drawnow;
            waitforbuttonpress
        end
        warning('on');

    end
%- STIM
    function plotMeanS(obj)
        fld=obj.getSFld;
        X=obj.(fld).X;
        S=obj.(fld).s;
        ctgInd=obj.(fld).ctgInd;
        nPix=size(S,1);
        n=numel(X);
        m=zeros(S,n,1);
        for i = 1:n
            m(:,i)=mean(S(:,ctgInd));
        end
        obj.fig('Mean Stim');
        for i = 1:n
            if i ==1
                hold off;
            end
            plot(m(:,i),'lineWidth',2);
            hold n;
        end
        axis square;
        Axis.format();
    end
%- FILTERS
    function plotF0(obj,fSmp,bNew)
        if nargin < 2
            fSmp=[];
        end
        if nargin < 3
            bNew=true;
        end

        f=obj.fig('AMA Filters',bNew,[1180 990]);

        bFFT=false;
        fld=obj.getOFld;
        AMA.plotFilters(obj.(fld).f0, fSmp, obj.nFset, obj.fSize, obj.nFSplit,bFFT);
        Fig.suptitle(['Filters ' newline obj.get_title()]);
    end
    function plotF(obj,fSmp,bNew)
        if nargin < 2
            fSmp=[];
        end
        if nargin < 3
            bNew=true;
        end

        f=obj.fig('AMA Filters',bNew,[1180 990]);

        bFFT=false;
        fld=obj.getOFld;
        AMA.plotFilters(obj.(fld).f, fSmp, obj.nFset, obj.fSize, obj.nFSplit,bFFT);
        Fig.suptitle(['Filters ' newline obj.get_title()]);
    end
    function plotFF(obj,fSmp,bNew)
        if nargin < 2
            fSmp=[];
        end
        if nargin < 3
            bNew=true;
        end

        f=obj.fig('AMA Filters',bNew);

        if bNew
            set(f,'Position',[f.Position(1:2) 1180 990]);
        end
        bFFT=true;
        AMA.plotFilters(obj.(fld).f, fSmp, obj.nFset, obj.fSize, obj.nFSplit,bFFT);
        Fig.suptitle(['FF Filters ' newline obj.get_title()]);
    end
    %- TSNE
        function plotTSNE(obj)
            obj.plot_TSNE('');
        end
        function plotTSNEs(obj)
            obj.plot_TSNE('s');
        end
        function fitTSNE(obj,fNs,distM)
            obj.fit_TSNE('');
        end
        function fitTSNEs(obj,fNs,distM)
            obj.fit_TSNE('s');
        end
        function replotTSNE(obj)
            obj.fit_TSNE('');
            obj.plot_TSNE('');
        end
        function replotTSNEs(obj)
            obj.fit_TSNE('s');
            obj.plot_TSNE('s');
        end
        %-TSNE HIDDEN
        end
        methods(Hidden)
            function fit_TSNE(obj,fld,fNs)
                if nargin < 3 || isempty(fNs)
                    fNs=1:obj.nF;
                elseif numel(fNs)==1
                    fNs=1:fNs;
                end
                if nargin < 3 || isempty(distM)
                    distM='cosine';
                end
                if strcmp(fld,'s')
                    [R,~]=obj.getResponses();
                else
                    [~,R]=obj.getResponses();
                end
                if strcmp(fld,'s')
                    fld2='tsnes';
                else
                    fld2='tsne';
                end
                obj.(obj.getOFld).(fld2)=tsne(R(:,fNs),'Distance',distM,'Standardize',true,'Perplexity',50);
            end
            function plot_TSNE(obj,fld)
                if strcmp(fld,'s')
                    fld2='tsnes';
                    nm=' simple';
                else
                    fld2='tsne';
                    nm='';
                end
                if ~isfield(obj.(obj.getOFld),fld2)
                    obj.fit_TSNE(fld);
                end

                f=obj.fig(['AMA Filters' nm]);

                hold off;
                Y=obj.(obj.getOFld).(fld2);
                hold off;
                gscatter(Y(:,1),Y(:,2),obj.S.ctgInd);
                Axis.format('R1','R2');
                axis square;
                Fig.suptitle(['t-SNE Responses ' newline obj.get_title()]);

            end
        end
        methods
%- WEIGHTS
    function plotW(obj,bNew)
        if nargin < 2
            bNew=true;
        end
        f=obj.fig('AMA Weights',bNew);
        imagesc(obj.OUT.w);
        Axis.format;
        axis square;
        Fig.suptitle(['AMA Weights ' newline obj.get_title()]);
        colormap hot;
        colorbar;
    end
    function plotWFF(obj,f)
        fld=obj.getOFld;
        if nargin < 2 || isempty(f)
            f=obj.(fld).f;
        end
        ff=RNorm.getff(f);
        W{1}=RNorm.frqW(f,ff);
        W{2}=RNorm.phsW(f,ff);
        W{3}=RNorm.cmplxW(f,ff);

        cl=Num.minMax([W{1}(:); W{2}(:); W{3}(:)]);
        t={'freq','phase','complx'};
        F=obj.fig('AMA FF Weights');
        for i = 1:length(W)
            subplot([1,3],i);
            hold off;
            imagesc(W{i});
            caxis(cl);
            axis square;
            Axis.format('','',t{i});
            colormap hot;
            colorbar;
        end
        Fig.sgtitle(['FF Weights ' newline obj.get_title()]);

    end
    function plotWDiff(obj,f)
        fld=obj.getOfld;
        if nargin < 2 || isempty(f)
            f=obj.(fld).f;
        end
        W{1}=RNorm.diffW(f);
        W{2}=RNorm.diffW2(f);
        W{3}=RNorm.diffW3(f);
        t={'diffW','diffW2','diffW3'};
        F=obj.fig('AMA diff Weights');
        for i = 1:length(W)
            subplot([1,3],i);
            hold off;
            imagesc(W{i});
            caxis(cl);
            axis square;
            Axis.format('','',t{i});
            colormap hot;
            colorbar;
        end
        Fig.sgtitle(['Diff Weights ' newline obj.get_title()]);
    end
    function plotWSim(obj,f)
        obj.getOFld;
        if nargin < 2 || isempty(f)
            f=obj.(fld).f;
        end
        W=RNorm.simWeights(f);
        WN=RNorm.simWeightsNeg(f);
        W2=RNorm.simWeights2(f);
        WN2=RNorm.simWeights2(f);

        cl=Num.minMax([W(:); WN(:); W2(:); WN2(:)]);

        f=obj.fig('AMA Sim Weights');
        subplot([2,2],1,1);
        hold off;
        imagesc(W);
        caxis(cl);
        axis square;
        Axis.format('','','W');
        colormap hot;
        colorbar;

        subplot([2,2],2, 1);
        hold off;
        imagesc(WN);
        caxis(cl);
        axis square;
        Axis.format('','','WN');
        colorbar;

        subplot([2,2],1, 2);
        hold off;
        imagesc(W2);
        caxis(cl);
        axis square;
        Axis.format('','','W2');
        colorbar;

        subplot([2,2],2, 2);
        hold off;
        imagesc(WN2);
        caxis(cl);
        axis square;
        Axis.format('','','WN2');
        colorbar;

        Fig.sgtitle(['Sim Weights ' newline obj.get_title()]);
    end
%- RESPONSES
    function plotRS(obj,fPairs,ctg2plt,axisLims)
        fld=obj.getOFld;
        if nargin < 2
            if size(obj.(fld).f,2)>1
                fPairs=[1 2];
            else
                fPairs=[1 1];
            end
        end
        if nargin < 3
            ctg2plt=[];
        end
        if nargin < 4
            axisLims=[];
        end
        plotT=0.2;

        f=obj.fig('AMA R Simple');
        hold off;
        r=obj.getResponses();
        AMA.plotResponses(r,obj.S.ctgInd,obj.S.X,fPairs,ctg2plt,true,false,axisLims,plotT);
        Fig.sgtitle(['Responses simple ' newline obj.get_title()]);
    end
    function plotREachCtg(obj,normType,pauseT,bNew)
        if nargin < 2
            normType=obj.normType;
        end
        if nargin < 3 || isempty(pauseT);
            pauseT=0.3;
        end
        if nargin < 4 || isempty(bNew);
            bNew=false;
        end
        REC=obj.getMREachCtg(normType);

        fld=obj.getOFld;
        fld2=obj.getSFld;
        X=obj.(fld2).X;

        f=obj.fig('R Each Ctg',bNew);
        hold off;
        yl=Num.minMax([real(REC(:)); imag(REC(:))]);

        for i = 1:size(REC,1)
            set(groot,'CurrentFigure',f.Number);
            plot(X,REC(i,:),'lineWidth',2);
            if i==1
                Axis.format('Ctg','R',['R-Norm ' normType]);
                ylim(yl);
            end
            hold on;

            if pauseT > 0
                drawnow
                pause(pauseT);
            end
        end
        strs=cellfun(@(x) ['f' num2str(x)],num2cell(1:obj.nF),'UniformOutput',false);
        legend(strs);
    end
    function plotRSumEachCtg(obj,normType)
        if nargin < 2
            normType=obj.normType;
        end
        REC=obj.getMREachCtg(normType);

        fld=obj.getOFld;
        fld2=obj.getSFld;
        X=obj.(fld2).X;

        f=obj.fig('R Sum Each Ctg');
        hold off;
        plot(X,mean(REC,1),'lineWidth',2);
        xlim(Num.minMax(X));
        Axis.format('Ctg','R Mean',['R-Norm ' normType]);
    end
    function plotRDistEach(obj,normType)
        if nargin < 2 || isempty(normType)
            normType=[];
        end
        [~,r]=obj.getResponses(normType);

        fld=obj.getOFld;
        fld2=obj.getSFld;
        X=obj.(fld2).X;
        ctgInd=obj.(fld2).ctgInd;
        n=numel(X);

        nF=size(r,2);
        RC=ceil(sqrt(nF));
        obj.fig('plotRDistEach');
        for i = 1:nF
            subplot(RC,RC,i);
            hold off;
            for j = 1:n
                [c,x]=hist(r(j==ctgInd,i));
                plot(x,c,'lineWidth',1);
                hold on;
                axis square;
                Axis.format('R mag','count',['f' num2str(i)]);
            end
        end

    end
    function plotRDist(obj,normType)
        if nargin < 2 || isempty(normType)
            normType=[];
        end
        [~,r]=obj.getResponses(normType);

        fld=obj.getOFld;
        fld2=obj.getSFld;
        X=obj.(fld2).X;
        ctgInd=obj.(fld2).ctgInd;
        n=numel(X);

        obj.fig('plotRDist');
        nF=size(r,2);
        for i = 1:nF
            if i == 1
                [c,x]=hist(r(:,i),50);
                hold off;
            else
                [c]=hist(r(:,i),x);
            end
            plot(x,c,'lineWidth',2);
            hold on;
            axis square;
            Axis.format('R mag','count',['f' num2str(i)]);
        end

    end
    function plotR(obj,normType,fPairs,ctg2plt,axisLims)
        fld=obj.getOFld;
        if nargin < 2 || isempty(normType)
            normType=obj.normType;
        end

        if nargin < 3 || isempty(fPairs);
            if size(obj.(fld).f,2)>1
                fPairs=[1 2];
            else
                fPairs=[1 1];
            end
        end
        if nargin < 4
            ctg2plt=[];
        end
        if nargin < 5
            axisLims=[];
        end

        plotT=0.2;
        ctgIndTst=[];

        f=obj.fig('AMA R Complex');
        hold off;
        [~,r]=obj.getResponses(normType);
        AMA.plotResponses(r, obj.S.ctgInd, obj.S.X,fPairs,ctg2plt,true,false,axisLims,plotT);
        Fig.sgtitle(['Responses ' newline obj.get_title()]);
    end
    function plotRFits(obj,fPairs,ctg2plt,axisLims)
        if nargin < 2
            fPairs=[];
        end
        if nargin < 3
            ctg2plt=[];
        end
        if nargin < 4
            axisLims=[];
        end

        % TODO?
        sTst=[];
        ctgIndTst=[];
        fTst=[];

        f=obj.fig('AMA R Spinner');
        hold off;
        AMA.plotResponses(r,obj.S.ctgInd,obj.X,[],sTst,fTst,ctgIndTst, ctg2plt,false,true,axisLims);
        Fig.sgtitle(['R Fits ' newline obj.get_title()]);
    end
%- ERRORS
    function plotErrorCumuCtg(obj,errorType,normType)
    end
    function plotErrorCumu(obj,errorType,normType)
        if nargin < 2
            errorType=[];
        end
        if nargin < 3
            normType=[];
        end
        Em=zeros(obj.nF,1);
        for i = 1:obj.nF
            Em(i)=obj.getError(errorType,normType,1:i);
        end

        obj.fig('ErrorCumu');

        plot(Em,'lineWidth',2);
        Axis.format('Filter No.','Error');
        axis square;
    end
    function plotErrorEachCtg(obj,eType,normType)
        if nargin < 2
            eType=[];
        end
        if nargin < 3
            normType=[];
        end
        if nargin < 4
            pauseT=[];
        end
        EEC=obj.getErrorEachCtg(eType,normType);

        fld=obj.getSFld;
        IN=obj.(fld);

        f=obj.fig('Error Each Ctg');
        hold off;
        set(groot,'CurrentFigure',f.Number);
        plot(IN.X,EEC,'lineWidth',2);
        hold on;
        Axis.format('Ctg','Error');
        %strs=cellfun(@(x) ['f' num2str(x)],num2cell(1:obj.nF),'UniformOutput',false);
        %legend(strs);
    end
    function plotErrorEachF(obj,eType,pauseT)
        if nargin < 2
            eType=[];
        end
        if nargin < 3
            pauseT=[];
        end
        EEF=obj.getErrorEachF(eType);

        fld=obj.getSFld;
        IN=obj.(fld);

        f=obj.fig('Error Each Ctg and Filter');
        hold off;
        for i = 1:obj.nF
            set(groot,'CurrentFigure',f.Number);
            plot(IN.X,EEF(i,:),'lineWidth',2);
            hold on;
            if pauseT > 0
                drawnow
                pause(.3);
            end
        end
        Axis.format('Ctg','Error');
        strs=cellfun(@(x) ['f' num2str(x)],num2cell(1:obj.nF),'UniformOutput',false);
        legend(strs);
    end
    function plotError(obj,errorType)
        if nargin < 2 || isempty(errorType)
            errorType=obj.errorType;
        end
        obj.getError(errorType);
        [counts,edges]=hist(obj.OUT.errorAll);

        f=obj.fig('AMA Errors');
        plot(edges,counts,'k','LineWidth',2);
        Axis.format('Error','Count');
        g=gca;
        set(g,'YScale','log');
        Axis.yticksLog(g);
        Fig.suptitle(['AMA Error ' errorType newline obj.get_title()]);
    end
%- VIEW
    function view(obj,order,inds)
        if nargin < 2 || isempty(order)
            order='worst';
        end
        switch order
        case 'worst'
            o='descend';
        case 'best'
            o='ascend';
        end
        vOpts=[];
        if isfield(obj.S,'ptchOpts') && ~isempty(obj.S.ptchOpts)
            [LX,RX]=obj.get_Pctrs();
        end
        [~,err]=obj.getError();
        XHat=obj.getEstiamtes();
        if nargin < 3 || isempty(inds)
            [~,inds]=sort(err,o);
        elseif ismember(0,inds);
            inds=find(inds);
        end

        LRy=[-.5 .5];
        i=1;
        if isfield(obj.S,'Ixyz')
            N=2;
        else
            N=1;
        end
        while i <= length(inds)
            ii=inds(i);
            X=obj.S.X(obj.S.ctgInd(ii));
            if ~isempty(LX)
                vOpts.Dsp=X;
                Lx=[LX(ii) LX(ii)];
                Rx=[RX(ii) RX(ii)];
            end
            [LCor,R]=obj.splitS(ii);
            obj.fig('AMA vieww worst',i==1,[],true);
            subPlot([1 N],1,1);
            hold off;
            plot(LCor,'r'); hold on
            plot(R,'b');
            if ~isempty(vOpts)
                plot(Lx,LRy,'b');
                plot(Rx,LRy,'r');
            end
            titl=[ 'error ' num2str(err(ii)) newline ...
                   'X '     num2str(X) newline...
                   'XHat '  num2str(XHat(ii)) ...
                ];
            Axis.format('','',titl);
            axis square;

            if N>1
                subPlot([1 N],1,2);
                plot(squeeze(obj.S.Ixyz(:,:,3,ii)));
            end
            Axis.format('','','');
            axis square;

            i=plotFlipper(i);
        end
    end
end
methods(Static)
%- STATIC Filters
%- STATIC RESPONSES
end
end
