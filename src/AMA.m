classdef AMA < handle & AMA_train & AMA_other
%: f    [ nPix x nF ]
%: PCor [ nStim x 1]
%: PAll [ nStim x nCtg]
%: S    [ nPix  x nStm ]
%: XXX LS
%:
%: FORMATS
%: Trn
%:    f00 - fully random, no params, no C
%:    f0  - fully random, no C
%:    fIn - no C
%: Index
%:    F - packed
%:    f - (unpacked)
%: (Saved)
%:    f - Complex
properties
    alias
    trainORtest=''
    initTime
end
% XXX read only
properties
    hash
    fname
end
properties
    Stim
    Nrn
    Obj %objective

    Opts
end
properties(Hidden)
    bCaller

    STORE % RM?
    bStored=false

    last=struct()
    bSet=false

    AMADeps={}
    bHasFit=0;

    % TMP
    tmpOpts
    bTEST
end
properties(Hidden)
%properties(Access=protected)
    StimTmp
    StimTrn
    StimTst
end
methods(Static)
%- TESTS
    function test0()
        AMA('bTEST',1);
    end
    function test1()
        obj=AMA('bTEST',1,'sTrnAlias','DSP21_0');
        obj.init();
    end
    function test2()
        obj=AMA('bTEST',1,...
                'sTrnAlias','DSP21_0',...
                 'nF',2,...
                 'nFset',1 ...
        );
        obj.train();
    end
    function test3()
        obj=AMA('bTEST',0,...
                'sTrnAlias','DSP21_0',...
                 'nF',2,...
                 'nFset',1 ...
        );
        obj.train();
    end
    function test4()
        obj=AMA('bTEST',0,...
                'sTrnAlias','DSP21_0',...
                'nF',4,...
                'nFset',4, ...
                'bWindow',true, ...
                'bNormalizeS',true ...
        );
        %obj.train();
    end
    function testNrw()
        obj=AMA('bTEST',0,...
                'sTrnAlias','DSP21_0',...
                'nF',2,...
                'nFset',2, ...
                'bAnalytic',true, ...
                'fCmplx',1, ...
                'normType','nrw', ...
                'bFourierF',true...
        );
        obj.train();
    end
    function testNrw2()
        obj=AMA('bTEST',0,...
                'sTrnAlias','DSP21_0',...
                'initTime','2023-04-21-14:44',...
                'nF',2,...
                'nFset',1, ...
                ...
                'bWindow',true, ...
                'bNormalizeS',true, ...
                ...
                'MaxIter',100, ...
                'nRecFit',1, ...
                ...
                'bAnalytic',true, ...
                'fCmplx',3, ...
                'bCmplxParam',1, ...
                'bSplitParam',1, ...
                ...
                'bFourierF',true,...
                'normType','nrw'  ...
        );
        %obj.loadP();
        %obj.loadF();
        %obj.train();
        %obj.recurseAll();
    end

%- Load Static
    function loadAlias(obj)
        % TODO
    end
    function out=isAlias(in)
        % TODO
        out=false;
    end
    function S=loadParams(in)
        if nargin < 1
            error('Argument required')
        end
        dire=AMA.getDir(in);
        fname=[dire 'params.mat'];

        if ~Fil.exist(fname)
            error('File does not exist')
        end

        S=load(fname);
    end
    function f=loadFilters(in,fF)
        if nargin < 2
            fF=obj.getDefaultfF();
            %error('Argument required')
        end
        name=AMA.fF2name(fF);
        dire=AMA.getDir(in);
        fname=[dire name '.mat'];

        if ~Fil.exist(fname)
            error('File does not exist')
        end
        load(fname);
    end
end
methods
    function obj=AMA(varargin)
        if nargin < 1
            return
        end
        obj.constr(varargin{:});
    end
    function obj=test(obj)
        %obj.Obj.test(obj.Stim);
    end
    function constr(obj,varargin)
        obj.mkChildren();

        % PARSE
        if (nargin == 1) && ischar(varargin{1});
            opts=AMA.loadOpts(varargin{1});
            if nargin == 2
                opts{'stimOpts'}{'trainORtest'}=varargin{2};
            end
            obj.parse(opts);
        else
            obj.parse(varargin);
        end

        Obj.assignin(obj,'AMA');
        obj.init();
    end
%- CONSTRUCT
    function mkChildren(obj)
        chlds=obj.ls_chlds();
        opts=obj.ls_opts();

        % populate children
        for i = 1:length(chlds)
            chld=chlds{i};
            if startsWith(chld,'Stim')
                chldStr='Stim';
            else
                chldStr=chld;
            end
            str=['AMA' chldStr '();'];
            %str
            obj.(chld)=eval(str);
        end

        % populate opts
        obj.Opts=struct();
        for i = 1:length(opts)
            opt=opts{i};
            obj.Opts.(opt)=eval(['AMA' opt 'Opts();']);
        end

        % send objects
        for i = 1:length(chlds)
            chld=chlds{i};
            obj.(chld).ama_receive(obj);
            %obj.(chld).ama_init();
        end
        for i = 1:length(opts)
            opt=opts{i};
            obj.Opts.(opt).ama_receive(obj);
            %obj.(chld).ama_init();
        end
    end
    function setAliasAsDate(obj)
        obj.alias=obj.initTime;
    end
    function setInitTime(obj)
        obj.initTime=Date.timeFilStr();
    end
    function stim_mk_clone(obj)
        % TRAIN
        obj.StimTrn=obj.StimTmp;
        obj.StimTst=obj.StimTmp.copy();
        obj.StimTmp=[];

        obj.StimTrn.trainORtest='train';
        obj.StimTrn.sTstAlias='';
        obj.StimTrn.sTstFName='';

        obj.StimTst.trainORtest='test';
        obj.StimTst.sTrnAlias='';
        obj.StimTst.sTrnFName='';
    end
%- PARAMS
    function parse(obj,vargs)
        if iscell(vargs)
            vargs=struct(vargs{:});
        elseif ischar(vargs) && AMA.isAlias(vargs)
            obj.loadAlias();
            return
        end

        [~,vargs]=Args.applyIf(obj,vargs);
        vargs=Args.struct2row(vargs);

        % OPTS
        opts=obj.ls_opts();
        for i = 1:length(opts)
            opt=opts{i,1};
            P=obj.Opts.(opt).getP();
            if isempty(P)
                continue
            end
            [O,vargs]=Args.simpleLoose([],P,vargs{:});
            vargs=Args.bicell2row(vargs);
            obj.Opts.(opt).ama_apply_opts(O);
        end

        % Children
        chlds=obj.ls_chlds();
        for i = 1:length(chlds)
            chld=chlds{i,1};
            P=obj.(chld).getP();
            if isempty(P)
                continue
            end
            [O,vargs]=Args.simpleLoose([],P,vargs{:});
            vargs=Args.bicell2row(vargs);
            obj.(chld).ama_apply_opts(O);
        end

        % DUPLICATE
        obj.stim_mk_clone();

    end
%- INIT
    function init(obj)
        obj.Stim.init();
        obj.Nrn.init();
        obj.Obj.init();

        obj.assert();
        if isempty(obj.initTime)
            obj.setInitTime();
        end
    end
    function assert(obj)
        typs=obj.ls_opts();
        for i = 1:length(typs)
            typ=typs{i};
            obj.Opts.(typ).assert();
        end
    end
    function P=getPs(obj,bVal,bSave,trainORtest)
        if nargin< 2 || isempty(bVal)
            bVal=false;
        end
        if nargin< 3 || isempty(bSave)
            bSave=false;
        end
        if nargin< 3 || isempty(trainORtest)
            trainORtest=obj.trainORtest;
        end
        %if nargin< 3
        %    SkipFlds={};;
        %end

        opts=obj.ls_opts();
        root=obj.Opts;
        P1=obj.getP_fun_(root,opts,bVal,bSave,trainORtest);

        opts=obj.ls_chlds();
        root=obj;
        P2=obj.getP_fun_(root,opts,bVal,bSave,trainORtest);

        P=[P1; P2];
    end
    function P=getP_fun_(obj,root,opts,bVal,bSave,trainORtest)
        P=[];
        for i = 1:length(opts)
            opt=opts{i,1};
            if strcmp(opt,'StimTmp')
                if strcmp(trainORtest,'train')
                    opt='StimTrn';
                elseif strcmp(trainORtest,'test')
                    opt='StimTst';
                elseif strcmp(trainORtest,'none')
                    continue
                end
                fopt='Stim';
            else
                fopt=opt;
            end
            %if ismember(opt,SkipFlds)
            %    continue
            %end

            if bSave
                p=root.(opt).getPSave();
            else
                p=root.(opt).getP();
            end
            if isempty(p)
                continue
            end
            if bVal
                flds=p(:,1);
                V=[];
                for j = 1:length(flds)
                    v={root.(opt).(flds{j})};
                    V=[V; v];
                end
                if bSave && startsWith(opt,'Stim')
                    p=obj.(opt).done2do();
                end
            else
                V=[];
            end
            f=repmat({fopt},size(p,1),1);
            pp=[f p V];
            pp(cellfun(@isempty,V),:)=[];
            P=[P; pp];
        end
    end
%- load
    function load(obj)
        obj.loadF();
        obj.loadP(obj);
    end
    function loadF(obj)
        obj.loadF_(2);
        obj.loadF_(-4);
    end
    %- Filters
    function loadF_(obj,fF)
        if nargin < 2
            fF=obj.getDefaultfF();
        end
        f=AMA.loadFilters(obj,fF);
        switch fF
        case {1,2}
            sz=size(f,2);

            obj.Nrn.f=f;
            obj.Nrn.nF=sz;

            obj.Trn.f=f;
            if isempty(obj.Opts.Trn.nF) || obj.Opts.Trn.nF < sz
                obj.Opts.Trn.nF=sz;
            end
        case -4
            obj.Trn.loadInF(f);
        case -1
            obj.Nrn.f0=f0;
        otherwise
            error('invalid fF')
        end
    end
    %- Params
    function loadP(obj,in)
        if nargin < 2 || isempty(in)
            in=obj;
        end
        P=AMA.loadParams(in);
        obj.constr(P);
    end
    %- Stim
    function loadS(obj,fnameORalias,trainORtest)
        if nargin < 1
            fnameORalias=[];
        end
        if nargin < 2
            trainORtest=[];
        end
        obj.Stim.load(fnameORalias,trainORtest);
    end
%- Save
    function save(obj)
        ef=obj.check_overwrite();
        if ef; return; end

        obj.saveF_(2,false); % Trn (OUT)
        obj.saveF_(-4,false); % FIn

        obj.saveParams_('train',false);
    end
    function saveF(obj)
        obj.saveF_(2,true); % Trn (OUT)
        obj.saveF_(-4,false); % FIn
    end
    function saveTrainParams(obj)
        obj.saveParams_('train',false);
    end
    %- check
    function exitflag=check_overwrite(obj);
        dire=AMA.getDir(obj);
        exitflag=true;
        if ~Dir.exist(dire) || isempty(Dir.files(dire))
            exitflag=false;
            return
        end
        out=Input.yn('Files exist. Overwrite?');
        if out==1
            exitflag=false;
            Dir.mk_p(dire);
        else
            exitflag=true;
        end
    end
    %-  filters
    function save_training_(obj)
    end
    function saveF_(obj,fF,bCheck)
        if nargin < 3 || isempty(bCheck)
            bCheck=true;
        end
        [f,fF]=obj.getF(fF);
        if bCheck
            exitflag=obj.check_overwrite();
            if exitflag; return; end
        end

        dire=AMA.getDir(obj);

        name=AMA.fF2name(fF);
        fname=[dire name];
        save(fname,'f');
    end
    %- params
    function saveParams_(obj,trainORtest,bCheck)
        if nargin < 3 || isempty(bCheck)
            bCheck=true;
        end
        P=obj.getPs(true,true,trainORtest);
        P(:,3:4)=[];
        if bCheck
            exitflag=obj.check_overwrite();
            if exitflag; return; end
        end

        obj.params2Yaml_(P);
        obj.params2Mat_(P);
    end
    function params2Mat_(obj,P)
        S=cell2struct(P(:,3),P(:,2));
        dire=AMA.getDir(obj);
        fname=[dire 'params'];
        save(fname,'-struct','S');
    end
    function params2Yaml_(obj,P)
        % TO STRUCT
        flds1=unique(P(:,1));
        S=struct();
        for i = 1:length(flds1)
            fld1=flds1{i};
            S.(fld1)=struct();

            inds=ismember(P(:,1),fld1);
            flds=P(inds,2);
            vals=P(inds,3);
            S.(fld1)=cell2struct(vals,flds);
        end

        dire=AMA.getDir(obj);

        fnamey=[dire '_params.yaml'];
        yaml.WriteYaml(fnamey,S,0);
        %TO YAML
        %XXX
    end
%- SET
    function unsetCaller(obj)
        obj.bCaller=false;
    end
    function out=get.Stim(obj,val)
        %cl=onCleanup(@() obj.unsetCaller());
        %obj.bCaller=true;
        trainORtest=obj.trainORtest;
        if isempty(trainORtest)
            bTrn=obj.StimTrn.bLoaded;
            bTst=obj.StimTst.bLoaded;
            if ~bTrn && ~bTst
                bTrn=~isempty(obj.StimTrn.fname);
                bTst=~isempty(obj.StimTst.fname);
            end
            if ~bTrn && ~bTst
                bTrn=~isempty(obj.StimTrn.alias);
                bTst=~isempty(obj.StimTst.alias);
            end
            if bTrn & bTst
                error('trainORtest needs to be set')
            elseif bTst
                trainORtest='test';
            elseif bTrn
                trainORtest='train';
            end
        end
        switch trainORtest
        case 'train'
            out=obj.StimTrn;
        case 'test'
            out=obj.StimTest;
        end
    end
    %- ALIAS
    function set.alias(obj,val)
        % XXX
    end
    %- CHANGE
    function set.trainORtest(obj,val)
        if ~any(strcmp(val,{'train','test'}))
            error('invalid trainORtest value')
        end
        if obj.bCaller
            obj.trainORtest=val;
        else
            obj.change_fun(val,false);
        end
        obj.bSet=true;
    end
    function reset(obj)
        % XXX
        % Reset main parameters to fit name

        if ~obj.bSet
            error('nothing to reset')
        end
        obj.alias=obj.last.alias; % XXX
        obj.last=struct();
    end
    function change(obj,trainORtest,bForce)
        if nargin < 3 || isempty(bForce)
            bForce=false;
        end
        if nargin < 2
            if strcmp(obj.trainORtest,'test')
                trainORtest='train';
            elseif strcmp(obj.trainORtest,'train')
                trainORtest='test';
            else
                error('invalid')
            end
            disp(sprintf('Changing to %s',trainORtest));
        end
        if ~any(strcmp(trainORtest,{'train','test'}))
            error('Invalid')
        end
        obj.change_fun(trainORtest);
    end
    function change_fun(obj,trainORtest,bForce)
        if nargin < 3
            bForce=false;
        end
        obj.bCaller=true;
        cl=onCleanup(@() cl_fun(obj));

        bTrainCur=strcmp(obj.trainORtest,'train');
        bTrainNew=strcmp(trainORtest,'train');

        % LOAD IF NEEDED
        if ~bTrainNew % TEST
            if ~obj.StimTst.bLoaded
                obj.StimTst.load();
            end
        elseif bTrainNew
            if ~obj.StimTrn.bLoaded
                obj.StimTrn.load();
            end
        end
        obj.trainORtest=trainORtest;

        function cl_fun(obj)
            obj.bCaller=false;
        end
    end
%- UTIL
    function [R,Nrn]=getR(obj,fF,fR,bCombine)
        if nargin < 2
            fF=[];
        end
        if nargin < 3 || isempty(fR)
            %fR='r';
            %fR='R';
            fR='Rm';
        end
        if nargin < 4 || isempty(bCombine)
            bCombine=false;
        end

        Nrn=obj.Nrn.copy();
        if nargin < 7 || isempty(f)
            f=obj.getF(fF);
        end
        Nrn.f=f;

        % GET RESPONSE
        Nrn.respond;
        R=Nrn.(fR);
        if bCombine
            nSplit=obj.Stim.nSplit;
            nF=Nrn.nF;
            N=size(R,1);
            % XXX
            p=[1 3 2];
            %p=[1 2 3];
            R=permute(reshape(R,N,nSplit,nF),p);
            R=sum(R,3);
        end


    end
    function fF=getDefaultfF(obj)
        if ~isempty(obj.Trn.f)
            % DEFAULT TO TRN.f IF TRAINED
            fF=2;
        else
            % DEFAULT TO Nrn.f IF TRAINED
            fF=1;
        end
    end
    function f=trnF2f_(obj,fld,bPack)
        if nargin < 3
            bPack=true;
        end
        T=obj.Trn.copy();
        T.fIn=T.(['f' fld]);
        T.cIn=T.(['c' fld]);
        f=complex(T.fIn,T.cIn);
        if ~bPack
            return
        end

        P=obj.Prg.copy();
        P.incIMax();

        I=obj.Index.copy();
        I.pack(T,P);

        %I.modifyBounds();
        %I.plotBounds();
        %[fraw,w,c]=obj.separate(fIn);

        f=I.unpack(I.F);
    end
    function [f,fF]=getF(obj,fF)
        if nargin < 2 || isempty(fF)
            fF=obj.getDefaultfF();
        end
        switch fF
        case -4
            % fIn
            f=obj.trnF2f_('In',true);
        case -3
            % f00
            f=obj.trnF2f_('00');
        case -2
            % NRN.f0
            f=obj.Nrn.f0;
        case -1
            % TRN.f0
            f=obj.trnF2f_('0');
        case 0
            % Index.f0
            f=obj.Index.F;
            if obj.Nrn.fCmplx > 0;
                f=complex(f(:,:,1),f(:,:,2));
            end
        case 1
            % NRN F
            f=obj.Nrn.f;
        case 2
            % Trn F
            f=obj.Trn.f;
        end
    end
%- Helpers
end
methods(Access=?AMAChld)
    function out=ls_opts(obj);
        out={...
             'Obj';
             'Trn';
             'FMin';
        };
    end
    function out=ls_chlds(obj)
        out={'StimTmp';
             'Nrn';
             'Obj';
             'Trn';
             'Prg';
             'Index';
             'Optim';
             'Plot';
        };
    end
    function ls_init_order()
        % XXX
        out={...
             'Obj';
             'Trn';
             'FMin';
        };
        out={'Stim';
             'Nrn';
             'Obj';
             'Trn';
             'Prg';
             'Index';
             'Optim';
             'Plot';
        };
    end
    function ls_copy(obj)
    end
end
methods(Static)
%- UTIL
    function out=fF2name(fF)
        switch fF
        case -4
            out='fIn';
        case -3
            out='f00';
        case -2
            %NRN
            out='f0';
        case -1
            out='f0';
        case 0
            out='F';
        case 1
            % NRN
            out='f';
        case 2
            % Trn
            out='f';
        end

    end
    function dire=getDirBase()
        dire=[Env.var('DATA') 'AMA' filesep];
    end
    function out=getDir(in)
        if isa(in,'AMA')
            if isempty(in.alias)
                alias=in.initTime;
            else
                alias=in.alias;
            end
        else
            alias=in;
        end
        out=[AMA.getDirBase() alias filesep];
    end
end
end
