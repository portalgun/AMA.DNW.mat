classdef AMAFMinOpts < handle & AMAChld & AMAOpts
properties

    Algorithm
    LargeScale
    SpecifyConstraintGradient
    FiniteDifferenceType
    ScaleProblem
    Display
    MaxFunEvals
    MaxIter
    TolFun
    TolX
    TolCon
    UseParallel

    btchSz
    stpSzMax
    stpSzMin
    stpSzEta

    bGlobal
    MaxTime
    GDisplay

    GradObj
    GradConstr
    OutputFcn

    BasinRadiusFactor
    DistanceThresholdFactor
    PenaltyThresholdFactor
    NumStageOnePoints
    NumTrialPoints

    %PlotFcn
end
properties(Access=protected)
    Optim
    Obj
    TrnOpts
    ObjOpts
end
properties(Hidden)
    AMADeps={'Optim','Obj','ObjOpts','TrnOpts'};
end
methods(Static)
    function P=getP();
        P=[AMAFMinOpts.getP_fmin();...
           AMAFMinOpts.getP_SGD();...
           AMAFMinOpts.getP_global();...
        ];
    end
    function P=getP_fmin()
        %alg='active-set';
        %alg='interior-point'; % sparse & large
        %alg='trust-region-reflective';
        alg='sqp'; % medium
        P={...
           'Algorithm',   alg,'';
           'LargeScale'   'off','';
           'SpecifyConstraintGradient', false,'';
           'FiniteDifferenceType', 'central',''; % default 'forward'
           'ScaleProblem', true,'';
           'Display',     'iter','';
           'MaxFunEvals', [],'';
           'MaxIter',     [],'';
           'TolFun',      1e-8,'';
           'TolX',        1e-8,'';
           'TolCon',      1e-6,'';
           'UseParallel', [],'';
        };
    end
    function P=getP_SGD()
        P={ ...
           'btchSz',      [],'';
           'stpSzMax',    [],'';
           'stpSzMin',    [],'';
           'stpSzEta',    [],'';
        };
    end
    function P=getP_global()
        P={ ...
           'bGlobal',                 [],'isBinary_e';
           'MaxTime',                 inf,''; %XXX
           'GDisplay',                'iter','';
            ...
           'BasinRadiusFactor',       0.2,'';
           'DistanceThresholdFactor', 0.75,'';
           'PenaltyThresholdFactor',  0.2,'';
           ...
           'NumStageOnePoints',       50,'';
           'NumTrialPoints',          250,'';
       };
    end
    function out=getPAliases()
        out={...
            'GDisplay','Display';
        };
    end
    function out=listFMinOpts();
        P=AMAFMinOpts.getP_fmin();
        out=P(:,1);
    end
    function out=listGlobalOpts()
        P=AMAFMinOpts.getP_global();
        out=P(:,1);
    end
    function out=listSGDOpts()
        P=obj.getP_SGD();
        out=P(:,1);
    end
end
methods
    function obj=AMAFMinOpts(varargin)
        if nargin < 1
            return
        end
    end
    function out=get.Optim(obj)
        out=obj.AMA.Optim;
    end
    function out=get.ObjOpts(obj)
        out=obj.AMA.Opts.Obj;
    end
    function out=get.TrnOpts(obj)
        out=obj.AMA.Opts.Trn;
    end
    function parse(obj,varargin)
    end
    function assert(obj)
        O=obj.ObjOpts;
        if any(strcmp(O.ppAlg,{'FLLGPU','GRDGPU'})) && ~isempty(obj.UseParallel) && obj.UseParallel
            disp(['Cannot use parallel with FLLGPU or GRDGPU']);
        end
    end
%- List
    function out=retFMin(obj)
        selflds=obj.listFMinOpts;
        out=obj.ret_fun_(selflds,[]);
    end
    function out=retGlobal(obj)
        selflds=obj.listGlobalOpts();
        out=obj.ret_fun_(selflds,[]);
    end
    function out=ret_fun_(obj,selflds,rmflds)
        obj.init();
        obj.assert();
        prps=obj.prps;

        prps=obj.prps;
        if ~isempty(selflds)
            prps=prps(ismember(prps,selflds));
        elseif ~isempty(rmflds)
            prps(ismember(prps,rmflds))=[];
        else
            error('selflds or rmflds are emtpy')
        end


        aliases=obj.getPAliases();
        out=struct();
        for i = 1:length(prps)
            prp=prps{i};
            if ismember(prp,aliases(:,1))
                ind=find(ismember(aliases(:,1),prp));
                aprp=aliases{ind,2};
            else
                aprp=prp;
            end
            out.(aprp)=obj.(prp);
        end
    end
%- INIT
    function opts=init(obj)
        T=obj.TrnOpts;


        if isempty(obj.bGlobal)
            obj.bGlobal=false; % TODO
        end

        if isempty(obj.MaxFunEvals) && ~obj.bGlobal
            %obj.MaxFunEvals=obj.nPix*800;
            obj.MaxFunEvals=10000 * max([1,T.nFset-1 * 10]) * max([1,T.nWset-1 * 10]);
        end

        if isempty(obj.MaxIter) && ~obj.bGlobal
            obj.MaxIter=100 * max([1,T.nFset-1 * 10]) * max([1,T.nWset-1 * 10]);
        end

        obj.OutputFcn=obj.Optim.getOutputFun();
        O=obj.ObjOpts;
        if any(strcmp(O.ppAlg,{'GRD','GRDGPU'}))
            obj.GradObj='on';
            obj.GradConstr='on';
        end
        if any(strcmp(O.ppAlg,{'FLLGPU','GRDGPU'}))
            obj.UseParallel = false;
        elseif isempty(obj.UseParallel)
            obj.UseParallel = true;
        end

    end
end
end
