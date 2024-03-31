classdef AMA_plot_models < handle
properties
end
methods(Static)
end
methods
    function plotModelBare(obj)
        obj.plot_model([],'none');
    end
    function plotModelBareM(obj)
        obj.plot_model([],'none',true);
    end
    function plotModel(obj,normType)
        if nargin < 2
            normType=strtok(obj.normType,'__');
        end
        obj.plot_model([],normType);
    end
    function plotModelM(obj,normType)
        if nargin < 2
            normType=strtok(obj.normType,'__');
        end
        obj.plot_model([],normType,true);
    end
    function plotModelAtten(obj,normType)
        if nargin < 2
            normType=strtok(obj.normType,'__');
        else
            normType=strtok(normType,'__');
        end
        normType=[normType '__atten'];
        obj.plot_model([],normType);
    end
    function plotModelHeeger(obj,normType)
        if nargin < 2
            normType=strtok(obj.normType,'__');
        else
            normType=strtok(normType,'__');
        end
        %normType=[normType '__broad'];
        obj.plot_model([],'heeger');
    end
    function plotModelBroad(obj,normType)
        if nargin < 2
            normType=strtok(obj.normType,'__');
        else
            normType=strtok(normType,'__');
        end
        %normType=[normType '__broad'];
        obj.plot_model([],'broad');
    end
    function plotModelNarrow(obj,normType)
        if nargin < 2
            normType=strtok(obj.normType,'__');
        else
            normType=strtok(normType,'__');
        end
        %normType=[normType '__broad'];
        obj.plot_model([],'narrow');
    end
    function plotModelBroadM(obj,normType)
        if nargin < 2
            normType=strtok(obj.normType,'__');
        else
            normType=strtok(normType,'__');
        end
        %normType=[normType '__broad'];
        obj.plot_model([],'broad',true);
    end
    function plotModels(obj,other)
        models={'none','broad','narrow','heeger','std','W1','multi'};
        if nargin >= 2 && ~isempty(other)
            models=[models other];
        end
        obj.plot_model([],models,false);
    end
%- MAIN
    function plot_model(obj,baseMod,normType,bMods)
        cl=onCleanup(@() cl_fun(obj));
        obj.bParentPlot=true;
        warning('off');
        if nargin < 4
            bMods=false;
        end

        %Nrms={'none','broad',normType};
        if ~iscell(normType)
            Nrms={normType};
        else
            Nrms=normType;
        end
        bMultiNrm=numel(Nrms)>1 & ~bMods;

        if bMods
            Mods={'none','abs','e','log','Log','2S','R2S','R2S_2S'};
            WH=[1085 670];
        elseif bMultiNrm
            bMultiNrm=true;
            Mods={'none'};
            WH=[1260 1175];
            if numel(normType)==5
                WH=[1260 1175];
            end
        else
            Mods={'none'};
            WH=[1200 290];
        end
        if ~isempty(baseMod)
            Mods=strcat(baseMod,'_',Mods);
        end

        F=numel(Nrms);
        R=5;
        C=numel(Mods);
        errors=zeros(F,C);
        ncorrect=zeros(F,C);
        maxRs=zeros(F,C);
        mRs=zeros(F,C);
        names=cell(F,C);
        for f = 1:F
        for m = 1:C

            nrm=Nrms{f};
            modu=Mods{m};
            name=[nrm '__' modu];
            names{f,m}=name;
            errors(f,m)=obj.getError([],name);
            [~,~,ncorrect(f,m)]=obj.getNCorrect('ME',name);

            resp=obj.getMREachCtg(name);
            mRs(f,m)=max(abs(resp(:)));

            [~,resp]=obj.getResponses(name);
            maxRs(f,m)=max(resp(:));
        end
        end
        bFlip=~bMods & ~bMultiNrm;
        if bMultiNrm
            C=F;
            F=1;
            names=names';
            errors=errors';
            ncorrect=ncorrect';
            maxRs=maxRs';
            mRs=mRs';
        end
        RC=[R C];

        mx=max(errors(:));
        yle=[0 mx*1.05];
        for f = 1:F

            Name=[strrep(obj.alias,'_','-') newline 'NormType ' strrep(Nrms{f},'_','-')];
            ff=Fig.new(Name);
            set(ff,'Position',[ff.Position(1:2) WH]);
            for m = 1:C
                c=0;
                name=names{f,m};
                if bMultiNrm
                    titl=strrep(name,'__none','');
                else
                    titl=strrep(Mods{m},'_',' ');
                end
                if m==1
                    titl2=sprintf('MaxR = %.3f       ',maxRs(f,m));
                else
                    titl2=sprintf('%.3f',maxRs(f,m));
                end

                %- R Dist
                c=c+1;
                subPlot(RC,c,m,bFlip);
                obj.plotRDist(name);
                xlabel('');
                title(titl);
                if m~=1
                    ylabel('');
                else
                    xlabel('R');
                end
                g=gca;
                set(g.Legend,'visible','off');

                %- R EACH CTG
                c=c+1;
                subPlot(RC,c,m,bFlip);
                obj.plotREachCtg(name,0);
                title(titl2);
                xlabel('');
                %title(strrep(Mods{m},'_',' '));
                if m~=1
                    ylabel('');
                else
                    ylabel('Mean R');
                    xlabel('Ctg');
                end
                ylim([-mRs(f,m) mRs(f,m)]);
                g=gca;
                set(g.Legend,'visible','off');

                %- R EACH CTG2
                c=c+1;
                subPlot(RC,c,m,bFlip);
                obj.plotREachCtg([name '_^2'],0);
                title('');
                xlabel('');
                %title(strrep(Mods{m},'_',' '));
                if m~=1
                    ylabel('');
                else
                    ylabel('Mean R^2');
                    xlabel('Ctg');
                end
                yl=max(abs(ylim));
                ylim([0 yl]);
                g=gca;
                set(g.Legend,'visible','off');


                %- RSUM
                c=c+1;
                subPlot(RC,c,m,bFlip);
                obj.plotRSumEachCtg(name);
                xlabel('');
                title('');
                if m~=1
                    ylabel('');
                else
                    ylabel('Mean R');
                    xlabel('Ctg');
                end
                yl=max(abs(ylim));
                ylim([-yl yl]);


                %- ERROR
                c=c+1;
                subPlot(RC,c,m,bFlip);
                obj.plotErrorEachCtg([],name);
                title('');
                if m~=1
                    ylabel('');
                end
                if m == 1
                    str1='Em =';
                    str2='C =';
                else
                    str1='';
                    str2='';
                end
                str1=sprintf('%s %.4f',str1,errors(f,m));
                str2=sprintf('%s %.4f',str2,ncorrect(f,m));
                xlabel(str1);
                title(str2);
                ylim(abs(yle));

                %xlabel('');
                %if i==2; ylabel(''); end
                %if i==2; ylabel(''); end
                %if i > 1
                %    xlabel('');
                %end
            end
            Fig.supTitle(Name);
        end
        %disp('Errors:')
        %for f = 1:F
        %for m = 1:C
        %    name=names{f,m};
        %    spl=strsplit(name,'__');
        %    base=spl{1};
        %    mods=strrep(spl{2},'_',' ');
        %    disp(sprintf('  %8s %6s --- %.3f',base,mods,errors(f,m)));
        %end
        %end

        function cl_fun(obj)
            obj.bParentPlot=false;
            warning('on');
        end
    end
%- CMP
    function modelCmp(obj,obj2,fNum)
        if nargin < 3
            fNum=[];
        end
        flds={'bNoise','alpha','rMax','TRN'};
        %vals={false,0,5.7};
        vals={false,0,5.7};
        for i = 1:length(flds)
            fld=flds{i};
            S1.(fld)=obj.(fld);
            S2.(fld)=obj2.(fld);
        end
        cl=onCleanup(@() cl_fun(obj,obj2,S1,S2));
        for i =1:length(flds)
            fld=flds{i};
            if strcmp(fld,'TRN')
                continue
            end
            obj.(fld)=vals{i};
            obj2.(fld)=vals{i};
        end
        if ~isempty(fNum)
            if fNum==1
                obj2.TRN=obj.TRN;
            elseif fNum==2
                obj.TRN=obj2.TRN;
            end
        end

        % ERRORS
        err=abs(errMES2./errMES);
        err1=~bC1;
        err2=~bC2;

        % LAGS

        % LumStats

        mdl2=(X2.MaxDiffLag);
        mdl1=(X1.MaxDiffLag);
        gdInds=bC1 & bC2;
        bdInds=~gdInds & (errMAP2>1 | abs(errMES2)>1);
        otInds=~gdInds & ~bdInds;
        bdInds1=~bC1 & ~gdInds;
        bdInds2=~bC1 & ~gdInds;

        %figure(99)
        %hist(X2.CRnk(1,:))
        %figure(100)
        %hist(X2.mxInds);
        %size(X2.lag)
        %hist(X2.mxInds)

        X=obj.S.X(obj.S.ctgInd);
        C1=X1.C;
        C2=X2.C;
        bX1=true;
        if bX1
            X0=X1;
            err0=errMES;
            mdl0=mdl1;
            bC0=bC1;
            C0=C1;
        else
            X0=X2;
            err0=errMES2;
            mdl0=mdl2;
            bC0=bC2;
            C0=C2;
        end


        figure(95)
        hold off;
        plot3(X(~bC2),lgs2(~bC2),errMES2(~bC2),'k.'); hold on;
        plot3(X(bC2), lgs2(bC2), errMES2(bC2),'b.'); hold on;
        plot3([0,-15],[0,8],[0 0],'LineWidth',3);
        Axis.format('X','lag rank error');
        %xlim([-10 10])
        zlabel('err');
        % only between mdl 10

        figure(96)
        hold off;
        LCor=Lu1.ml-Lu1.mr;
        plot3(X(~bC1),LCor(~bC1),errMES(~bC1),'k.'); hold on;
        plot3(X(bC1), LCor(bC1), errMES(bC1),'b.'); hold on;
        %plot3([0,-15],[0,8],[0 0],'LineWidth',3);
        Axis.format('X','lum');

        return
        while I <=size(C0,2)
            figure(95)
            i=sinds(I);

            %-
            Elag=X1.ELag(i,1);
            err=errMES(i,1);
            XHat=err+X(i);
            Olag=X1.MaxOLag(i);

            subPlot([1 2],1,1);
            hold off;
            r=X1.lagRnk(:,i);
            y=(C1(X1.CRnk(:,i),i));

            plot(r,'.k');
            hold on;
            plot([0,63],[Elag,Elag],'r');
            plot([0,63],[Olag,Olag],'b');

            Axis.format(['Max LCor ' num2str(Olag) newline ...
                         'Exp LCor ' num2str(Elag)  newline...
                         'X '     num2str(X(i))          newline ...
                         'XHat  ' num2str(XHat)          newline ...
                         'Error ' num2str(err)  ...
            ]);
            title('X1');

            %-
            Elag=X2.ELag(i,1);
            err=errMES2(i,1);
            XHat=err+X(i);
            Olag=X2.MaxOLag(i);

            subPlot([1 2],1,2);
            hold off;
            r=X2.lagRnk(:,i);
            y=(C2(X2.CRnk(:,i),i));

            plot(r,'.k');
            hold on;
            plot([0,63],[Elag,Elag],'r');
            plot([0,63],[Olag,Olag],'b');
            Axis.format(['Max LCor ' num2str(Olag) newline ...
                         'Exp LCor ' num2str(Elag)  newline...
                         'X     ' num2str(X(i))          newline ...
                         'XHat  ' num2str(XHat)          newline ...
                         'Error ' num2str(err)  ...
            ]);
            title('X2');

            I=plotFlipper(I);
        end
        return

        return

        figure(96)
        hold off;
        plot(X2.MaxOLag,errMES,'k.'); hold on;

        %- GD BAD OTHER
        figure(98)
        hold off;
        plot(errMAP(bdInds),errMAP2(bdInds),'.r');hold on
        plot(errMAP(gdInds),errMAP2(gdInds),'.b');
        plot(errMAP(otInds),errMAP2(otInds),'.k');

        %- errMAP sep
        figure(99)
        hold off;
        inds1=bdInds1 & errMAP>2;
        inds2=bdInds2 & errMAP2> 2.4; %errMAP2>3.3369;
        plot(errMAP(bdInds),errMAP2(bdInds),'.k'); hold on
        plot(errMAP(inds1),errMAP2(inds1),'r.');
        plot(errMAP(inds2),errMAP2(inds2),'m.');

        %plot(mdl1(bdInds),mdl1(bdInds1),'r.'); hold on;
        %plot(mdl1(inds3),mdl2(bdInds2),'k.'); hold on;
        return


        %- max lag 1
        figure(100)
        bdUp=bdInds & ~inds2 & ~inds1;
        hold off;
        %plot3(mdl2(ginds),errMES(ginds),errMES2(ginds),'k.');
        inds3=bdInds1 & abs(mdl1)>1;
        plot(mdl1(bdInds),mdl2(bdInds),'r.'); hold on;
        plot(mdl1(inds3),mdl2(inds3),'k.'); hold on;

        %- max lag 2
        figure(101)
        hold off;
        bdUp=bdUp & ~inds3;
        %plot3(mdl2(ginds),errMES(ginds),errMES2(ginds),'k.');
        inds4=bdInds2 & abs(mdl2)>1;
        plot(mdl1(bdInds),mdl2(bdInds),'r.'); hold on;
        plot(mdl1(inds4),mdl2(inds4),'k.'); hold on;

        % lumdiff
        figure(102)
        hold off;
        bdUp=bdUp & ~inds4;

        inds5=bdUp & bdInds2;
        plot(errMAP(inds5),abs(X2.MeanOCorr(inds5)),'.');
        %plot(Lu2.vB(bdInds2 & bdUp),errMAP2(bdInds2 & bdUp),'.k');

        %plot3(ludiff(bdInds),errMES(bdInds),errMES2(bdInds),'k.');
        %hold on
        %plot3(ludiff(bdUp),errMAP(bdUp),errMAP2(bdUp),'b.');

        %- LAST
        figure(900)
        hold off;
        plot(errMES,errMES2,'k.'); hold on;
        plot(errMES(inds1),errMES2(inds1),'r.');
        plot(errMES(inds2),errMES2(inds2),'r.');
        plot(errMES(inds3),errMES2(inds3),'r.');

        figure(901)
        lstInds=bdInds & ~inds1 & ~inds2 & ~inds3 & ~inds4;
        plot(errMES(lstInds),errMES2(lstInds),'k.');

        %plot3(mdl2(inds),errMES(inds),errMES2(inds),'b.');
        %plot3(mdl2(inds2),errMES(inds2),errMES2(inds2),'r.');
        %plot3(mdl2(inds3),errMES(inds3),errMES2(inds3),'m.');
        %plot3(mdl2(inds4),errMES(inds4),errMES2(inds4),'g.');
        %inds3=abs(mdl1) > 1 & abs(mdl2) < 2 & errMES > 1.5;
        %inds4= gdInds & abs(errMES) > 1;

        return

        %plot(errMES(~inds),errMES2(~inds),'b.');
        %plot(errMES(inds),errMES2(inds),'r.');
        if false
            figure(99)
            hold off
            plot3(mdl1(~inds3),mdl2(~inds3),errMES(~inds3),'k.');
            hold on
            plot3(mdl1(inds3),mdl2(inds3),errMES(inds3),'r.');
            plot3(mdl1(inds4),mdl2(inds3),errMES(inds4),'r.');
            Axis.format('lag diff 1','lag diff 2');
            zlabel('errMES1');
        end



        Axis.format('Dif Max Lag','err MSE');
        zlabel('errMES2');

        hold on;
        Axis.format('Bi lum diff','err MES');
        zlabel('errMES2');

        figure(102)
        hold off;
        %errs=sort(err);
        %ind2=errs>=1;
        %ind1=errs<1;
        %ind0=errs==1;
        %e=sum(ind1);
        %e2=sum(ind2);
        hist(err,50);
        %plot(1:e,errs(ind1),'k.');
        %hold on;
        %plot(e+1:e+e2,errs(ind2),'r.');
        set(gca,'XScale','log');

        %xlim([-6 6]);
        %ylim([-6 6]);

        function cl_fun(obj1,obj2,S1,S2)
            Flds=fieldnames(S1);
            for j = 1:length(Flds)
                Fld=Flds{j};
                obj1.(Fld)=S1.(Fld);
                obj2.(Fld)=S2.(Fld);
            end
        end
    end
    function modelCmp2(obj1,obj2)
        obj1.fig('Bi cmp');
        obj1.bNoise=false;
        obj2.bNoise=false;
        %obj1.getStats();
        %obj2.getStats();
        %obj1.errAnalysis();
        %obj2.errAnalysis();
        %obj1.getBiDiff();
        %obj2.getBiDiff();

        obj2.TRN.f=obj1.TRN.f;
        obj2.getStats;

        obj1.fig('MES Cmp');
        hold off;
        plot(obj2.A.Em.X,obj2.A.Em.MES ,'.');
        hold on;
        plot(obj1.A.Em.X,obj1.A.Em.MES ,'.');

        i1=obj2.A.inds1;
        i2=obj2.A.inds2;
        i0=~i1 & ~i2;
        X=obj2.A.Em.X;

        obj1.fig('MES Cmp 4');
        Em=obj1.A.Em.MES;
        X=X-mean(X);
        %plot(Em+X/2,'.')

        obj1.fig('MES Cmp 3');
        hold off;
        %plot(X(i0),abs(obj1.A.Em.MES(i0)-obj2.A.Em.MES(i0)),'.');
        %plot(X(i0),(obj1.A.Em.MES(i0)-obj2.A.Em.MES(i0)),'.');
        %hold on;
        Ed=obj1.A.Em.MES-obj2.A.Em.MES;
        Em=obj2.A.Em.MES;
        Bid=obj1.A.bi-obj2.A.bi;
        Md=mean(obj1.A.P2.MSC-obj2.A.P2.MSC,2);
        ld=obj1.A.LCor.diffM-obj2.A.LCor.diffM;
        %plot3(X(i0),Ed(i0),Bid(i0),'k.'); hold on
        dsp2=obj2.A.dsp2;
        %ii=obj1.S.ctgInd==1;
        ii=true(obj1.S.nStim,1);
        corr(Md(ii & i0),Ed(ii & i0)) % .1692
        corr(Md(ii & i1),Ed(ii & i1)) % . 3370
        %for i = 1:obj.S.X
        %    plot3(X(i0),Bid(i0),Ed(i0),'k.'); hold on
        %    plot3(X(i1),Bid(i1),Ed(i1),'r.');
        %    zlabel('Bi diff');
        %    zlabel('Bi diff');
        %end

        obj1.fig('MES Cmp 2');
        hold off;
        plot(obj1.A.bi(i0)-obj2.A.bi(i0),abs(obj1.A.Em.MES(i0)-obj2.A.Em.MES(i0)),'.');
        hold on;
        plot(obj1.A.bi(i2)-obj2.A.bi(i2),abs(obj1.A.Em.MES(i2)-obj2.A.Em.MES(i2)),'.');
        %xlim([0 .25]);
        %ylim([0 .25]);


        %xlim([-.1,.1])
        %xlim([-.1,.1])
        %plot3(obj2.A.Em.X, obj2.A.dsp2, ,'.');
    end
end
end
