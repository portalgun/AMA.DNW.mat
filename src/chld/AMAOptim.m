classdef AMAOptim < handle & AMAChld
properties
    Name

    last
    best
    bestC

    state=''
    optimValues=struct('fval',inf,'constrviolation',inf)
    x % NOTE should be an optim value
    history % x history

    ind
    fInds
    wInds
    I

    bFMin
end
properties(Access=protected)
    TrnOpts
    Trn
    Prg
end
properties(Hidden)
    prevName=''
    bUpdate=true

    AMADeps={'TrnOpts','Trn','Prg'}
end
methods(Static)
end
methods
%- GET
    function out=get.TrnOpts(obj)
        out=obj.AMA.Opts.Trn;
    end
    function out=get.Trn(obj)
        out=obj.AMA.Trn;
    end
    function out=get.Prg(obj)
        out=obj.AMA.Prg;
    end
%-
    function obj=AMAOptim(varargin)
        if nargin < 1 || isempty(name)
            name='main';
        end
        obj.Name=name;
    end
    function init(obj)
        if bUpdate & strcmp(obj.Name,'main');
            obj.create_children();
            obj.bUpdate=false;
        end
    end
    function create_children(obj)
        flds={'last','best','bestC'};
        for i = 1:length(flds)
            obj.(flds{i})=AMAOptim([flds{i}]);
        end
        if obj.bContinue();
            obj.prompt_continue();
        end
        obj.history=[];
    end
    function bOut=bContinue(obj)
        bOut=obj.TrnOpts.bRecurse && strcmp(obj.state,'interrupt');
    end
    function swap(obj,cname)
        if isempty(cname) || strcmp(cname,'main')
            return
        end
        flds=obj.getFlds();
        for i = 1:length(flds)
            fld=flds{i};

            pTmpVal=obj.(fld);
            cTmpVal=obj.(cname).(fld);

            obj.(fld)=cTmpVal;
            obj.(cname).(fld)=pTmpVal;
        end
        obj.setState('sel',cname);
    end
%- util
    function cnames=getCNames(obj)
        cnames={'last','best','bestC'};
        cnames(~isfield(obj,flds))=[];
    end
    function [flds,inds]=getFlds(obj)
        flds={'state', ...
              'x','optimValues',...
              'history','fInds','wInds','I'};
        inds=[1,-1,-1,0,0,0,0];
    end
%- pack
    function [fw,fval,eflg,outpt,lambda,grd]=unpack(obj,fld)
        if isempty(fld) || strcmp(fld,'main')
            optm=obj;
        else
            optm=obj.(fld);
        end
        fw=optm.x;
        fval=optm.optimValues.fval;
        eflag=[]; % TODO construct?
        outpt=[]; % TODO from where? perhaps redundanct
        lambda=optm.optimValues.labmda;
        grd=optm.optimValues.gradient;
    end
    function setState(obj,state,fld)
        switch state
        case 'sel'
            obj.state=['selected ' fld];
            obj.prevName=fld;
        end
    end
%- user interrrupt
    function start(obj)
        obj.fInds=obj.Prg.fainds;
        obj.wInds=obj.Prg.wainds;
        obj.I=obj.Prg.getI();
        obj.bFMin=true;
    end
    function stop(obj)
        obj.bFMin=false;
    end
    function clfun(obj)
        if obj.bFMin % user iterrupt
             % ???
        end
    end
%- interrupt
    function clear(obj)
        [flds,inds]=obj.getFlds();
        for i =1:length(flds) % skipping state, x, optimValues, history
            if inds(i)==-1
                continue
            elseif inds(i)==1
                obj.(flds{i})='';
            else
                obj.(flds{i})=[];
            end
        end
        if strcmp(obj.Name,'main')
            flds={'last','best','bestC'};
            for i = 1:length(flds)
                obj.(fld{i}).clear();
            end
        else
            obj.clear_optim_values();
        end
    end
    function clearInterrupt(obj)
        obj.clearInterrupt();
    end
    function clear_optim_values()
        obj.x=[]; % NOTE X place inside optimValues?
        obj.optimValues=struct;
        obj.optimValues.fval=inf;
        obj.optimValues.constrviolation=inf;
    end
%- PROMPTS
    function prompt_continue(obj)
        if obj.TrnOpts.bSkipPromptInterrupt
            out=false;
        else
            out=Input.yn('Continue last interrupt?');
        end
        if out
            obj.saveInterrupt();
        else
            obj.clearInterrupt();
        end
    end
    function outFld=interruptPrompt(obj)
        cnames=obj.getCNames();
        n=length(cnames);

        if n==0
            outFld='';
            disp('Nothing to save')
            return
        elseif n==1
            outFld=cnames{1};
            disp('Single option to save')
            return
        end

        [f,c,x]=obj.collectOptimValues();

        % PRINT
        sstr='%d %s\n      fval   %.5f\n      constr %.8f';
        for i = 1:n
            str=sprintf(sstr, i,cnames,f(i),c(i));
            disp(str);
        end

        % INPUT
        out=Input.range(n);
        if out==0
            return
        end

        outFld=flds{out};
    end
    function collectOptimValues()
        obj.getCNames();
        n=length(cnames);

        z=zeros(n,1);
        f=z;
        c=z;
        x=cell(n,1);
        for i = 1:length(flds)
            F=obj.(flds{i});
            f(i)=F.optimValues.fval;
            c(i)=F.optimValues.constrviolation;
            x{i}=F.x;
        end
    end
%- FMINCON FUNS
    function fun=getOutputFun(obj)
        fun=@(x,optimValues,state) obj.output_fun(x,optimValues,state);
    end
    function fun=getPlotFun(obj)
        fun=@(x,optimValues,state) obj.plot_fun(x,optimValues,state);
    end
    function fun=getDispFun(obj)
        fun=@(x,optimValues,state) obj.disp_fun(x,optimValues,state);
    end
    function stop=output_fun(obj,x,optimValues,state);
        % SAVE STUFF FROM FMINCON

        % state init interrupt iter done
        stop=false;
        obj.state=state;

        % LAST
        obj.last.x=x;
        obj.last.optimValues=optimValues;

        % BEST FVAL
        if optimValues.fval < obj.best.optimValues.fval
            obj.best.optimValues=optimValues;
            obj.best.x=x;
        end

        % BEST CONSTR
        if optimValues.constrviolation < obj.bestC.optimValues.constrviolation
            obj.bestC.optimValues=optimValues;
            obj.bestC.x=x;
        end

        % HISTORY
        if isequal(state,'iter')
            if size(x,2) ~= size(obj.history,2);
                obj.history=x;
            else
                obj.history=[obj.history; x];
            end
        end
    end
    function stop=plot_fun(obj,x,optimvalues,state);
        % TODO
        stop=false;
        if mod(optimValues.funccount,10000)==0
            f=Fig.new('AMA Filters','bClear',false,'bFocus',false);
            AMA.plotFilters(obj.last.x, [], obj.TrnOpts.nFset, obj.TrnOpts.fSize, obj.TrnOpts.nFSplit);
            drawnow
        end
    end
    function disp_fun(obj,x,optimvalues,state)
        % TODO
        % Iter
        %       iteration
        % F-Count
        %       funccount
        % f(x)
        %       fval
        % Max constraint
        %       constrviolation
        % Line search steplength
        %       lssteplength
        % Directional derivative
        %       directionalderivative
        % First-order optimalityProceudre
        %       firstorderopt
        % Procedure
        %       procedure
    end
end
end
