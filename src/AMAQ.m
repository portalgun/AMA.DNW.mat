classdef AMAQ < handle
properties
    stAlias
    normType
    amaType
    fitType
    bMeanCon

    otherArgs

    baseName

    lst
    lstArgs

    errors
    prog
    lastMeth
end
properties(Constant)
    args={'stAlias','normType','fitType','bMeanCon','amaType'};
end
methods(Static)
    function out=script()
        base='DSP21_1_';
        dire='/Volumes/Data/.daveDB/AMA/';
        [names,fullnames]=Dir.files(dire);
        inds=startsWith(names,base) & endsWith(names,'.mat');
        files=fullnames(inds);
        names=strrep(names(inds),'.mat','');

        errors=inf(size(names));
        badInd=false(size(names));
        for i = 1:length(files)
            load(files{i});
            if ~isfield(ama,'OUT')
                %errors(i)=ama.TRN
                if isfield(ama.TRN,'f')
                    ind=find(~isempty(ama.TRN.error),1,'last');
                    errors(i)=ama.TRN.error(ind);
                else
                    badInd(i)=true;
                    continue
                end
            else
                errors(i)=ama.OUT.error(end);
            end
        end
        [errors,inds]=sort(errors);
        errors
        out=[names(inds) num2cell(errors)];
    end
    function out=getArg(args,name)
        out=args(ismember(AMAQ.args,'normType'));
    end
    function OUT=loadList();
        out=Cfg.readScript('AMAQueue');
        a=AMAQ.args;
        for i = 1:length(a)
            if ~isfield(out,a{i});
                out.(a{i})='';
            end
        end
        OUT=struct();
        OUT.baseArgs={out.stAlias,out.normType,out.fitType,out.bMeanCon,out.amaType};
        if isfield(out,'otherArgs')
            OUT.otherArgs=out.otherArgs;
        end
        OUT.lst=out.lst;
        OUT.lstArgs=out.lstArgs;
    end
end
methods
    function obj=AMAQ(baseArgs,lst,lstArgs,otherArgs)
        if nargin < 1
            return
        end
        obj.setBase(baseArgs{:});
        obj.lst=lst;
        obj.lstArgs=lstArgs;
        if nargin >= 4 && isempty(otherArgs)
            if iscell(otherArgs)
                obj.otherArgs=struct(otherArgs{:});
            elseif isstruct(otherArgs)
                obj.otherArgs=otherArgs;
            else
                error('Invalid otherArgs')
            end
        end
    end
    function setBase(obj,stAlias,normType,fitType,bMeanCon,amaType)
        if nargin == 2 && isstruct(stAlias)
            args=struct2cell(stAlias);
            obj.setBase(args{:});
            return
        elseif nargin ==2 && iscell(stAlias)
            obj.setBase(stAlias{:});
            return
        end
        % BASE PARAMS
        obj.stAlias=stAlias;
        obj.amaType=amaType;
        obj.fitType=fitType;
        obj.normType=normType;
        obj.bMeanCon=bMeanCon;
        obj.baseName=AMA.getAlias(obj.stAlias,obj.normType,obj.fitType,obj.bMeanCon,obj.amaType);
    end
    function Obj=setInst(obj,Obj,stAlias,normType,fitType,bMeanCon,amaType)
        Obj.stAlias=stAlias;
        Obj.fitType=fitType;
        Obj.normType=normType;
        Obj.bMeanCon=bMeanCon;
        Obj.amaType=amaType;
    end
    function out=loadBase(obj)
        [out,bSuccess]=AMA.load(obj.baseName);
        if ~bSuccess
            error('Base file does not exist')
        end
    end
    function out=loadInstExist(obj,args)
        name=AMA.getAlias(args{:});
        [out,bSuccess]=AMA.load(name);
        if ~bSuccess
            out=AMA.load(obj.baseName);
            out=obj.setInst(out,args{:});
        end
    end
    function out=getBaseArgs(obj)
        out={obj.stAlias,obj.normType,obj.fitType,obj.bMeanCon,obj.amaType};
    end
    function [args,name]=getInstArgs(obj,i)
        ind=ismember(AMAQ.args,obj.lstArgs);
        args=obj.getBaseArgs();
        [args{ind}]=obj.lst{i,:};
    end
    function Obj=setInstOtherArgs(obj,Obj,more)
        if ~isempty(obj.otherArgs)
            flds=fieldnames(obj.otherArgs);
            if numel(flds) > 0
                for i = 1:length(flds)
                    Obj.(flds{i})=obj.otherArgs.(flds{i});
                end
            end
        end
        if nargin < 3 || isempty(more) || numel(fieldnames(more)) < 1
            return
        end

        flds=fieldnames(more);
        for i = 1:length(flds)
            Obj.(flds{i})=more.(flds{i});
        end
    end
    function loadListFun(obj)
        OUT=AMAQ.loadList();
        obj.lst=OUT.lst;
        obj.lstArgs=OUT.lstArgs;
        if isfield(OUT,'baseArgs')
            obj.setBase(OUT.baseArgs);
        end
        if isfield(OUT,'otherArgs')
            obj.otherArgs=OUT.otherArgs;
        end
    end
    function record(obj,Obj,i)
        obj.errors(i)=Obj.OUT.error(end);
    end
%-- FIT
    function fitNew(obj,varargin)
        obj.main('new',varargin{:});
    end
    function fitRec(obj,varargin)
        obj.main('recurse',varargin{:});
    end
    function fitAll(obj,varargin)
        obj.main('recurseAll',varargin{:});
    end
    function fitFull(obj,varargin)
        obj.main('recurseFull',varargin{:});
    end
    function fitWeights(obj,varargin)
        obj.main('newW',varargin{:});
    end
    function cont(obj,varargin)
        if isempty(obj.lastMeth)
            error('Nothing to continue')
        end
        obj.main([],varargin{:});
    end
    function main(obj,meth,MaxIter)
        if ~isempty(meth)
            % NEW
            obj.lastMeth=meth;
            obj.errors=zeros(length(obj.lst),1);
            obj.prog=zeros(length(obj.lst),1);
            I=1:length(obj.lst);
            obj.loadListFun();
        else
            % CONTINUE
            meth=obj.lastMeth;
            I=find(obj.prog==1):length(obj.lst);
        end
        more=struct();
        if nargin >= 3
            more.fminconOpts.MaxIter=MaxIter;
        end
        bStart=true;
        for i = I
            obj.prog(i)=1;
            args=obj.getInstArgs(i);

            if strcmp(meth,'fitW')
                normType=AMA.getArg(args,'normType');
                if ~startsWith(normType,'W')
                    continue
                end
            end

            Obj=obj.loadInstExist(args);
            obj.setInstOtherArgs(Obj);

            assignin('base','obj',Obj);
            try
                Obj.(meth)();
            catch ME
                if bStart
                    rethrow(ME);
                end
                disp(ME.message);
                obj.prog(i)=-1;
                continue
            end
            obj.prog(i)=2;
            Obj.learnW();
            obj.record(Obj,i);
            bStart=false;
        end
    end
%----
    %function load(obj,name)
    %    for i = 1:obj.n
    %        name=[baseName '_' lst{i}];
    %        obj=AMA.load(name);
    %        disp(name);
    %        input('Press Return...');
    %        break
    %    end
    %end
    %function loadAll(obj)
    %    for i = 1:length(lst)
    %        name=[baseName '_' lst{i} '_all'];
    %        obj=AMA.load(name);
    %        disp(name);
    %        input('Press Return...');
    %        break
    %    end
    %end
    function error(obj)
        for i = 1:length(lst)
            name=[baseName '_' lst{i}];
            obj=AMA.load(name);
            disp(name);
            disp(obj.OUT.error(end));
            disp('');
        end
    end
    function errorBase(obj)
        for i = 1:length(lst)
            obj=AMA.load(baseName);
            obj.normType=lst{i};
            obj.getError;
            name=[baseName '_' lst{i}];
            disp(name);
            disp(obj.OUT.error(end));
            disp('');
            input('Press Return...');
        end
    end
    function plotF(obj)
        for i = 1:length(lst)
            name=[baseName '_' lst{i}];
            obj=AMA.load(name);
            obj.plotF;
            input('Press Return...');
        end
    end
    function plotR(obj)
        for i = 1:length(lst)
            name=[baseName '_' lst{i}];
            obj=AMA.load(name);
            obj.plotR;
            input('Press Return...');
        end
    end
end
end
