classdef AMA_train < handle
properties
    Trn
    Prg
    Index
    Optim
    Plot
end
properties(Hidden)
    Rst % TRN PROMPT RESET
    fRst=0
    bSuccessTrn
    lastStr
    gs
end
methods(Static)
end
methods
%- INIT
%- RST
    function cl=rst_init(obj)
        obj.bSuccessTrn=false;
        obj.fRst=1;
        obj.Rst=struct();
        obj.Rst.trainORtest=obj.trainORtest;
        obj.Rst.TrnOpts=obj.Opts.Trn.copy([],'AMATrnOpts');
        obj.Rst.Prg=obj.Prg;
        obj.lastStr='';
        cl=onCleanup(@() obj.rst_cleanup());
    end
    function rst_cleanup(obj)
        if ~obj.fRst || obj.bSuccessTrn
            return
        end
        obj.trainORtest=obj.Rst.trainORtest;
        obj.Opts.Trn=obj.Rst.TrnOpts;
        obj.Prg=obj.Rst.Prg;

        switch obj.fRst
        case 2
            obj.Trn=obj.Rst.Trn;
        end

        obj.fRst=0;
        obj.Rst=struct();
    end
    function rst_inc(obj)
        fR=obj.fRst+1;

        switch fR
        case 2
            obj.Rst.Trn=obj.Trn.copy();
        end

        obj.fRst=fR;

    end
%- Procedures
    function cont(obj,varargin)
        obj.train('bContinue',false,'bRecurse',false,'bAll',false,varargin{:});
    end
    function recurse(obj,varargin)
        obj.train('bContinue',false,'bRecurse',true,'bAll',false,varargin{:});
    end
    function recurseAll(obj,varargin)
        obj.train('bContinue',false,'bRecurse',true,'bAll',true,varargin{:});
    end
    function addF(obj,n,varargin)
        if nargin < 2 || isempty(n)
            n=1;
        end
        obj.Opts.Trn.nF=obj.Opts.Trn.nF+n;
        %obj.train('bRecurse',true,'bAll',true,varargin{:});
        obj.train('bContinue',true,'bRecurse',false,'bAll',false,varargin{:});
    end
%-
    function train(obj,varargin)
        % bAll = fix all but one at a time

        % wMode
        %   0 - no weight fitting
        %   1 - new weights
        %   2 - recurse
        %   3 - hyper
        %   4 - new hyper

        obj.init();
        if isempty(obj.trainORtest)
            obj.trainORtest='train';
        end

        % SAVE
        cl=obj.rst_init();

        % CHANGE
        obj.change('train');

        %% TRNOPTS
        obj.Trn.parseTrn(varargin{:});

        % PROMPT
        if obj.Prg.bPromptContinue
            % XXX
            %exitflag=obj.train_prompt();
            exitflag=false;
            if exitflag
                return % NOTE AUTO CLEANUP
            end
        end
        %% PRG
        obj.Prg.init_iter(); %BEFORE init_trn
        obj.fRst=2;

        %% TRN
        obj.rst_inc();
        obj.Trn.checkAndApply();

        % FIRST PLOT
        %obj.Plot.trn();

        %% Obj
        fun=obj.Obj.getFun();

        % CONS and BOUNDS
        [con,LB,UB]=obj.Trn.retCons();
        obj.Index.LB=LB;
        obj.Index.UB=UB;

        % fopts
        fOpts=obj.Opts.FMin.retFMin();
        gOpts=obj.Opts.FMin.retGlobal();

        obj.dispInit();
        while true
            [I,breakflag]=obj.Prg.incI();
            if breakflag; break; end

            % RECURSIVELY FIT SAME FILTERS nRecFit times
            while true
                [J,breakflag]=obj.Prg.incJ();
                if breakflag; break; end

                % INDECES AND F0
                obj.Index.pack(obj.Trn,obj.Prg);

                % BOUNDS
                [f0,lb,ub]=obj.Index.retF0();
                obj.Index.plotBounds();

                if J==1
                    obj.Plot.F0();
                    drawnow
                end

                % DISP
                obj.dispIter(I,J);

                % MINIMIZE
                exitflag=obj.minimize(fun,f0(:),lb(:),ub(:),con,fOpts,gOpts);

                % UPDATE
                obj.Prg.updateRecCompletion();

                % PLOT
                obj.Plot.FTrn();
                drawnow;

            end
            if breakflag; breakflag=false; end
            %if exitflag; break; end

            % UDPATE
            obj.Prg.updateCompletion();

            % SAVE
            if obj.bTEST
                continue
            end

            % SAVE
            if ~obj.bSuccessTrn % XXX HANDLE BETTER?
                obj.mark_success();
            end
            obj.save_training_();

            % PRG
            obj.Prg.disp_minTime();

        end
        obj.summary();
        obj.Prg.stop();
    end
    function [exitflag]=minimize(obj,fun,f0,lb,ub,con,fOpts,gOpts)
        OUT=cell(1,6);
        cl=onCleanup(@() obj.Optim.clfun());
        exitflag=true;

        A=[];
        b=[];
        Aeq=[];
        beq=[];

        problem=createOptimProblem('fmincon', ...
                                   'objective',fun, ...
                                   'x0',f0, ...
                                   'lb',lb, ...
                                   'ub',ub, ...
                                   'Aineq',A,...
                                   'bineq',b,...
                                   'Aeq',Aeq,...
                                   'beq',beq,...
                                   'nonlcon',con, ...
                                   'options',fOpts  ...
        );
        if gOpts.bGlobal && obj.Opts.Trn.bNewF
            obj.gs = GlobalSearch('options',gOpts);
        end

        obj.Optim.start();
        if obj.bTEST
            [OUT{:}]=obj.testMakeOut(fun,f0);
        elseif gOpts.bGlobal
            [OUT{1:5}]=run(obj.gs,problem);
        else
            [OUT{:}] = fmincon(problem);
        end
        exitflag=OUT{3};
        obj.Optim.stop();

        % REPACK
        obj.Trn.pack(obj.Prg.getI(),obj.Prg,OUT{:});
    end
    function [f,fval,exitflag,output,lambda,grad]=testMakeOut(obj,fun,f0);
        f=f0;
        fval=fun(f0);
        exitflag=false;
        output=struct(); % TODO
        lambda=nan; % TODO
        grad=nan; % TODO
    end
    function save_interrupt(obj,fld)
        % GET
        OUT=cell(1,6);
        [OUT{:}]=obj.optim.unpack(fld);

        % SET
        % % XXX recover other pack vals from PRG?
        obj.Trn.pack(I,[],OUT{:});

        % UPDATE OPTIM
        obj.optim.swap(fld);

        % XXX CHECK THESE
        obj.mark_success(true);
        obj.save_training_(); % SHOULD SAVE?
    end
    function dispInit(obj)
        ldiv=20;
        nSpc=2;

        div=repmat('-',1,ldiv);

        str0=sprintf('AMA TRAIN START%s\n',div);
        str1=obj.disp_params(nSpc);
        str2=obj.Prg.disp_info(nSpc);
        str=[str0 str1 newline str2];
        disp(str)
    end
    function out=disp_params(obj,nSpc)
        if nargin < 2 || isempty(nSpc)
            spc='';
        else
            spc=repmat(' ',1,nSpc);
        end

        % BASE MODEL
        str1=sprintf('%sData         %s\n',spc,obj.Stim.alias);
        str2=sprintf('%sNorm         %s\n',spc,obj.Nrn.normType);
        str3=sprintf('%sPosterior    %s\n',spc,obj.Opts.Obj.ppAlg);
        str4=sprintf('%sError        %s\n',spc,obj.Opts.Obj.errType);
        str5=sprintf('%sMean Con.    %d\n',spc,obj.Opts.Trn.bMeanCon);
        str=[str1 str2 str3 str4 str5];
        if nargout > 0
            out=str;
        else
            disp(str)
        end
    end
    function dispIter(obj,I,J)
        nSpc=2;
        I='';
        %I=append(I,obj.Prg.disp_train());
        I=append(I,obj.Obj.disp_error());
        I=append(I,obj.Obj.disp_nCorrect());
        I=append(I,obj.Nrn.disp_maxR());
        if obj.bTEST
            I=append(I,obj.Prg.disp_TEST(nSpc));
            obj.print(I);
        else
            obj.lastStr=I;
            obj.print(I);
        end
        function I=append(I,str)
            I=[I newline str];
        end
    end
    function print(obj,str)
        bDel=false;
        lastStr=obj.lastStr;
        bstr='';
        if ~isempty(lastStr) && bDel
            N=length(lastStr)-sum(lastStr=='\');
            bstr=repmat('\b',1,N);
            %fprintf(bstr);
        end
        str=strrep(str,newline,'\n');
        fprintf([bstr str]);
        obj.lastStr=str;
        if obj.bTEST
            out=input('','s');
            obj.lastStr=[obj.lastStr '\n' out];
        end
    end
    function mark_success(obj)
        obj.bSet=false;
        obj.bSuccessTrn=true;
        obj.bHasFit=1;
        obj.last=struct();
    end
%- INTERRUPT
    function saveInterrupt(obj)
        outFld=obj.optim.interruptPrompt();
        if isempty(outFld)
            return
        end
        obj.save_interrupt(outFld);
    end
    function saveInterruptLast(obj)
        obj.save_interrupt('last');
    end
    function saveInterruptBest(obj)
        obj.save_interrupt('best');
    end
    function saveInterruptBestConstr(obj)
        obj.save_interrupt('bestC');
    end
    function clearInterrupt(obj)
        obj.optim.clearInterrupt;
    end
%- PROMPT
    function testPrompt(obj)

        %disp(sprintf('TESTING: nIter %d of %d',i,obj.Prg.nIter));
        %obj.disp_minTime(i,obj.Prg.FLRN{i},obj.Prg.WLRN{i});
        input(sprintf('  PRESS RETURN'));
    end
    function exitflag=train_prompt(obj)
        % XXX
        exitflag=false;
        alias=obj.genAlias();
        if ~isempty(obj.alias) && obj.bHasFit && ~strcmp(obj.alias,alias)
            out=obj.prompt_change(alias);
            if ~out
                exitflag=true;
                return
            end
        end
        if obj.existSave(alias) ~isempty(obj.Trn) && ~isempty(obj.loadedAlias) && strcmp(obj.loadedAlias,alias)
            if (obj.T.bNewF) || (obj.T.bNewW)
                out=obj.prompt_exist();
                switch out
                case 'e'
                    disp('Exiting')
                    exitflag=true;
                    return
                case 'o'
                    return
                case 'l'
                end
            else
                return
                obj2=AMA.loadAlias();
                if isempty(obj2.complete)
                    return
                end
                out=obj.prompt_prg(obj2);
                switch out
                case 'e'
                    disp('Exiting')
                    exitflag=true;
                    return
                case 'r'
                    return
                case 'c'
                    obj.Prg=obj2.Prg;
                    obj.T=obj2.T;
                    obj.Trn=obj2.TRN;
                    obj.T.bContinue=true;
                end
            end
        end
    end
    function summary(obj)
    end
end
methods(Static)
%- RANDOM PARAMS
%- CONSTRAINTS
end
end
