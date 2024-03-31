classdef AMATrnOpts < handle & AMAChld & AMAOpts
properties
%- PARAMS
    bContinue
    bSkipInterruptPrompt
    bTEST
    bPlot

    nF
    nFset
    nWset
    FFix
    WFix
    nRecFit

    rndSd

    %- FILTER
    bLrnF
    bRecurse
    bAll
    bAllOne
    inds
    rinds
    excl
    nums

    %- FANCY
    bAnalytic
    bCmplxParam
    bSplitParam

    %- Weights
    bLrnW
    bRecurseW
    bAllW
    bAllOneW
    indsW
    rindsW
    exclW
    numsW


% other opts
    bMeanCon % XXX this doesn't belong here?

%- PARAMS COMPUTED
    bNewF
    bNewW

    bW

    fSz
end
properties(Access=protected)
    fCmplx
    ObjOpts
    Stim
    Nrn
end
properties(Hidden)
    bParsed=false
    rndSds % 1 f0, 2 w0, 3 other
    AMADeps={'ObjOpts','S','Nrn'};
end
methods(Static)
    function P=getP()
    %function train(obj,bRecurse,bAll,nums,excl,inds,rInds, wMode)
        constTest= @(x) isempty(x) || (numel(x)==1 & isnumeric(x));
        fTest    = @(x) isempty(x) || (numel(x)==1 && mod(x,1)==0 && x >  0);
        fTest0    = @(x) isempty(x) || (numel(x)==1 && mod(x,1)==0 && x >=  0);
        P={...
            'bMeanCon',        false, 'Num.isBinary'; % XXX MOVE?
            ...
            'bSkipInterruptPrompt',false,'';
            'bContinue',false,'';
            ...
            'bTEST',false,'';
            'bPlot',true,'';
            'nRecFit',  [],constTest;
            'nF',   [], fTest;
            'nFset',[],fTest;
            'nWset',0,fTest0;
            'FFix',[], fTest;
            'WFix',[], fTest;
            'rndSd',              [], constTest;
            ...
            'bLrnF',true,'';
            'bRecurse',false,'';
            'bAll',false,'';
            'bAllOne',false,'';
            'inds',[],'';
            'rinds',[],'';
            'excl',[],'';
            'nums',[],'';
            ...
            'bLrnW',false,'';
            'bRecurseW',false,'';
            'bAllW',false,'';
            'bAllOneW',false,'';
            'indsW',[],'';
            'rindsW',[],'';
            'exclW',[],'';
            ...
            'bAnalytic',0,'isBinary_e';
            'bCmplxParam',0,'isBinary_e';
            'bSplitParam',0,'isBinary_e';
        };
    end
end
methods
    function out=get.ObjOpts(obj)
        out=obj.AMA.Opts.Obj;
    end
    function out=get.Stim(obj)
        out=obj.AMA.Stim;
    end
    function out=get.fCmplx(obj)
        out=obj.AMA.Nrn.fCmplx;
    end
    function obj=init(obj)
        % nRecFit
        if isempty(obj.nRecFit)
            if strcmp(obj.ObjOpts,'SGD')
                obj.nRecFit=1;
            else
                obj.nRecFit=5;
            end
        end

        % master sd
        if isempty(obj.rndSd)
            obj.rndSd=randi(2^16,1);
        end

        % secondary seeds
        rng(obj.rndSd);
        obj.rndSds=randi(2^16,1,3);

        n=1;
        if obj.bLrnW
            n=n+1;
        end
        if obj.fCmplx > 0
            n=n+1;
        end
        if n==0
            n=[];
        end
        if isempty(obj.nF) && ~isempty(obj.AMA.Nrn.nF)
            obj.nF=obj.AMA.Nrn.nF;
        end
        if isempty(obj.nFset)
            obj.nFset=1;
        end
        if ~isempty(obj.nRecFit)
            obj.nRecFit=min([obj.nF obj.nRecFit]);
        end
        obj.fSz=[obj.Stim.nPix  obj.nF n]; % NOTE STIM
    end
    function parse(obj,varargin)
        Args.simple(obj,AMATrnOpts.getP,varargin{:});
        obj.bParsed=true;
    end
    function parseTrn(obj,varargin)
        Args.simpleExlusive(obj,AMATrnOpts.getP(),varargin{:});
    end
end
end
