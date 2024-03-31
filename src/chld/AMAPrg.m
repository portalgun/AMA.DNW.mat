classdef AMAPrg < handle & AMAChld
properties

    i
    I
    j
    iC
    jC
    nJ
    nI

    nIter
    ttic
    date

    complete
    completeW

    FLRN
    FFIX
    WLRN
    WFIX

    minTimeSec

    finds
    winds
    ffinds
    wfinds
    fainds
    wainds

    nIterF
    nIterW
    FI
    WI
end
% XXX SET READ ONLY
properties(Access=protected)
    nF
    nFset
    nWset
end
properties(Hidden)
    AMADeps={'TrnOpts'}
    Opts
end
% T
%
% bTest
% bNewF
% bNewW
% bContinue
%
% nRecFit
% nFset
% nWset
%
% inds
% excl
%
% con
% fOpts
%
% bPlot
%
methods
    function obj=AMAPrg(obj,T)
        if nargin < 1
            return
        end

    end
%- GET
    function out=get.Opts(obj)
        out=obj.AMA.Opts.Trn;
    end
    %
    function out=get.nF(obj)
        out=obj.Opts.nF;
    end
    function out=get.nFset(obj)
        out=obj.Opts.nFset;
    end
    function out=get.nWset(obj)
        out=obj.Opts.nWset;
    end
%- SET
    function set.nF(obj,nF)
        obj.Opts.nF=nF;
    end
    function set.nFset(obj,in)
        obj.Opts.nFset=in;
    end
    function set.nWset(obj,in)
        obj.Opts.nWset=in;
    end
%- Ret
    function I=getI(obj)
        if obj.i == -1
            I=-1;
        else
            I=obj.I(obj.i);
        end
    end
    function J=getJ(obj)
        J=obj.j;
    end
%- I
    function [I,breakflag]=incI(obj,val)
        if isempty(obj.i)
            if obj.nIter==0
                error('No iteration')
            end
            obj.iC=0;
            obj.nI=obj.nIter;
            obj.tic();
        elseif obj.i == -1
            % IF MORE HAVE BEEN ADDED
            obj.i=obj.iC;
            obj.nI=obj.nIter;
        end
        obj.j=[];
        breakflag=false;

        if nargin < 2 || isempty(val)
            i=obj.iC+1;
        else
            i=val;
        end
        if i > obj.nI
            i=-1;
            breakflag=true;
        end
        obj.i=i;
        I=obj.getI();

        if breakflag
            return
        end

        obj.finds=obj.FLRN{i};
        obj.winds=obj.WLRN{i};
        obj.ffinds=obj.FFIX{i};
        obj.wfinds=obj.WFIX{i};

        obj.fainds=unique([obj.finds obj.ffinds]);
        obj.wainds=unique([obj.winds obj.wfinds]);
    end
    function obj=updateCompletion(obj)
        obj.complete= [obj.complete  obj.FLRN{obj.i}];
        obj.completeW=[obj.completeW obj.WLRN{obj.i}];
        obj.iC=obj.i;
    end
%- J
    function [J,breakflag]=incJ(obj)
        if isempty(obj.j)
            obj.jC=0;
            obj.j=0;
            obj.nJ=obj.Opts.nRecFit;
        end
        breakflag=false;

        j=obj.jC+1;
        if j > obj.nJ
            j=-1;
            breakflag=true;
        end
        obj.j=j;
        J=obj.getJ();
    end
    function obj=updateRecCompletion(obj)
        obj.jC=obj.j;
    end
%- OTHER
    function incIMax(obj)
        obj.incI(obj.nI);
    end
%- Prompt
    function out=bPromptContinue(obj);
        out=~isempty(obj.i) && obj.i~=-1 && ...
            ~obj.Opts.bContinue && ~obj.Opts.bSkipInterruptPrompt;
    end
%- TIME
    function tic(obj)
        obj.date=date();
        obj.ttic=tic;
    end
    function t=toc(obj)
        t=toc(obj.ttic);
        obj.minTimeSec=t;
    end
    function stop(obj)

    end
%- DISP
    function out=disp_info(obj,nSpc)
        if nargin < 2 || isempty(nSpc)
            spc='';
        else
            spc=repmat(' ',1,nSpc);
        end

        % LEARN F INFO
        if ~isempty(obj.finds)
            n=Range.str(obj.finds);
        else
            n='-';
        end
        str1=sprintf('Learning Fs  %s\n',n);
        if ~isempty(obj.ffinds)
            n=Range.str(obj.ffinds);
        else
            n='-';
        end
        str2=sprintf('Fixing   Fs  %s\n',n);

        % LEARN W INFO
        if ~isempty(obj.winds)
            n=Range.str(obj.winds);
        else
            n='-';
        end
        str3=sprintf('Learning Ws  %s\n',n);
        if ~isempty(obj.wfinds)
            n=Range.str(obj.wfinds);
        else
            n='-';
        end
        str4=sprintf('Fixing Ws    %s\n',n);

        str=[spc str1 spc str2 spc str3 spc str4];

        str1=sprintf('nFset        %d\n',obj.nF);
        str2=sprintf('nRecFit      %d\n',obj.Opts.nRecFit);
        str=[str newline spc str1 spc str2];

        if nargout > 0
            out=str;
        else
            disp(str)
        end
    end
    function out=disp_minTime(obj)

        time=round(obj.minTimeSec);
        fstr=Range.str(obj.finds);
        wstr=Range.str(obj.winds);
        str1=sprintf('  Completed iter %d of %d in %d seconds:\n',i,obj.nIter,time);
        str2=sprintf('    Filters %s\n',fstr);
        str3=sprintf('    Weights %s\n',wstr);

        fstr2=Range.str(obj.complete);
        wstr2=Range.str(obj.completeW);
        str4=sprintf('  Total Progress\n');
        str5=sprintf('    Filters %s\n',fstr2);
        str6=sprintf('    Weights %s\n',wstr2);
        str=[str1 str2 str3 str4 str5 str6];
        if nargout > 0
            out=str;
        else
            disp(str);
        end
    end
    function out=disp_TEST(obj,nSpc)
        % XXX nRecFit
        if nargin < 2 || isempty(nSpc)
            spc='';
        else
            spc=repmat(' ',1,nSpc);
        end
        I='';
        I=append(I,sprintf('TESTING: nRec  %d of %d',obj.getJ,obj.Opts.nRecFit));
        I=append(I,sprintf('%sFLrn %s',spc,Range.str(obj.finds)));
        I=append(I,sprintf('%sFFix %s',spc,Range.str(obj.ffinds)));
        I=append(I,sprintf('%sWLrn %s',spc,Range.str(obj.winds)));
        I=append(I,sprintf('%sWFix %s',spc,Range.str(obj.wfinds)));
        I=append(I,sprintf('Press RETURN\n'));
        if nargout > 0
            out=I;
        else
            disp(I);
        end
       function I=append(I,str)
            I=[I newline str];
        end
    end
%- INIT
%
    function init_iter(obj)
        obj.Opts.init();
        obj.init_iterF();
        obj.init_iterW();

        % BOTH
        bRevF=obj.Opts.bAllW && ~obj.Opts.bAll  && ~obj.Opts.bAllOne;
        bRevW=obj.Opts.bAll  && ~obj.Opts.bAllW && ~obj.Opts.bAllOneW;

        u=unique([obj.WI obj.FI]);
        if isempty(u)
            error('Nothing to iterate (Hint: set nF or nW')
        end
        %obj.Prg.nIter=max([obj.Prg.nIterW obj.Prg.nIterF]);
        obj.nIter=numel(u);
        obj.I=u;

        obj.nIter=numel(obj.I);

        FFIX=cell(obj.nIter,1);
        FLRN=cell(obj.nIter,1);
        finds=obj.Opts.inds; % XXX
        nFset=obj.nFset;

        WFIX=cell(obj.nIter,1);
        WLRN=cell(obj.nIter,1);
        winds=obj.Opts.indsW;
        nWset=obj.nWset;

        if isempty(finds) && isempty(winds)
            error('empty finds and winds')
        end
        if nFset==0 && nWset==0
            error('Nothing to iterate (Hint: set nFset or nWset')
        end

        N=obj.nIter;
        for i = 1:N
            ffix=[];
            inds=(1:nFset)+(i-1)*nFset;
            inds(inds > obj.nF)=[];
            if ~isempty(finds)
                t=finds(inds);
            else
                t=[];
            end
            if obj.Opts.bAll
                FLRN{i}=finds;
            elseif ismember(i,obj.FI)
                FLRN{i}=t;
                if obj.Opts.bAllOne
                    ffix=finds(~ismember(finds,FLRN{i}));
                elseif bRevF
                    ffix=finds(~ismember(finds,finds(inds)));
                elseif min(inds)-1 > 0
                    ffix=finds(1:(min(inds)-1));
                end
            end
            FFIX{i}=unique([ffix obj.Opts.FFix]);
            FLRN{i}=FLRN{i}(~ismember(FLRN{i},FFIX{i}));

            wfix=[];
            inds=(1:nWset)+(i-1)*nWset;
            inds(inds > obj.nF)=[];
            w=winds(inds);
            if ~isempty(winds)
                w=winds(inds);
            else
                t=[];
            end
            if obj.Opts.bAllW
                WLRN{i}=w;
            elseif ismember(w,obj.WI)
                WLRN{i}=w;
                if obj.Opts.bAllOneW
                    wfix=winds(~ismember(winds,WLRN{i}));
                elseif bRevW
                    wfix=winds(~ismember(winds,winds(inds)));
                elseif min(inds)-1 > 0
                    wfix=winds(1:(min(inds)-1));
                end
            end
            WFIX{i}=unique([wfix obj.Opts.WFix]);
            WLRN{i}=WLRN{i}(~ismember(WLRN{i},WFIX{i}));
        end
        obj.FFIX=FFIX;
        obj.FLRN=FLRN;

        if ~isempty(obj.Opts.numsW)
            inds=obj.Opts.numsW;
        else
            inds=1:numel(WLRN);
        end
        obj.nIter=numel(inds);
        obj.WFIX=WFIX(inds);
        obj.WLRN=WLRN(inds);

        if obj.Opts.bContinue
            obj.i=obj.iC;
            obj.nI=obj.nIter;
        elseif obj.Opts.bRecurse
            obj.i=[];
            obj.j=[];
        end
    end
    function [p,inds,rInds]=init_inds(obj,i)
        % XXX OLD?

        %i=iternumber
        inds=1:obj.nFset+obj.nFset(i-1);

        %(i-1)*obj.Opts.nFset

        % FILTER INDICES TO BE LEARNED (NOT TO EXCEED nF)
        if ~isempty(obj.Opts.inds)
            p=1;
            inds=obj.Opts.inds;
        else
            p = obj.nFfix + obj.Prg.p;
            p(p > obj.nF)=[];
            inds   = (p+1):p+obj.nFset;
            inds(inds > obj.nF)=[];
        end

        % FILTERS TO BE FIXED

        if ~isempty(obj.Opts.rInds)
            rInds=obj.Opts.rInds;
        elseif obj.Opts.bAllOne
            rInds=1:obj.nF;
            rInds([obj.Opts.inds obj.Opts.excl])=[];
        else
            rInds=1:p;
        end

    end
    function init_iterF(obj)
        INDS=1:obj.Opts.nF;

        % F INDS
        %if isempty(obj.Opts.inds)
            obj.Opts.inds=INDS;
        %end

        %if ~isempty(obj.Opts.excl) || ~isempty(obj.Opts.FFix)
            obj.Opts.inds(ismember(obj.Opts.inds,[obj.Opts.excl obj.Opts.FFix]))=[];
        %end
        nF=numel(obj.Opts.inds);

        if obj.Opts.bAll
            obj.nFset=nF;
        elseif isempty(obj.nFset)
            obj.nFset=0;
        end

%F NUMS
        if obj.nFset==0
            obj.nIterF=0;
        else
            obj.nIterF=ceil(obj.nF/obj.nFset);
        end
        obj.Opts.bNewF=false;
        if ~obj.Opts.bRecurse && ~obj.Opts.bContinue
            obj.Opts.bNewF=true;
            obj.complete=[];
        elseif isempty(obj.complete) || ~obj.Opts.bContinue;
            obj.complete=[];
        end
        s=min(obj.Opts.inds(~ismember(obj.Opts.inds,obj.complete))); % XXX

        obj.FI=s:obj.nIterF;
    end
    function init_iterW(obj)
        INDS=1:obj.Opts.nF;

        % W INDS
        %if isempty(obj.Opts.indsW)
            obj.Opts.indsW=INDS;
        %end
        %if ~isempty(obj.Opts.exclW)
            obj.Opts.indsW(ismember(obj.Opts.indsW,[obj.Opts.exclW obj.Opts.WFix]))=[];
        %end
        nF=numel(obj.Opts.indsW);

        if obj.Opts.bAllW
            obj.nWset=nF;
        elseif isempty(obj.Opts.nWset)
            obj.nWset=1;
        end

        % W NUMS
        if obj.nWset==0
            obj.nIterW=0;
        else
            obj.nIterW=ceil(nF/obj.nWset);
        end
        obj.Opts.bNewW=false;

        if ~obj.Opts.bRecurseW && ~obj.Opts.bContinue
            obj.Opts.bNewW=true;
            obj.completeW=[];
        elseif isempty(obj.completeW) || ~obj.Opts.bContinue;
            obj.completeW=[];
        end
        s=min(obj.Opts.indsW(~ismember(obj.Opts.indsW,obj.completeW))); % XXX

        obj.WI=s:max(s,obj.nIterW);

    end
end
end
