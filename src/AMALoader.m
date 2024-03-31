classdef AMALoader < handle
properties
end
methods(Static)
    function fname=toOutFName(alias,dt)
        if nargin < 2 || isempty(dt)
            dt=Date.timeFileStr();
        end
        dire=[Env.var('DATA') 'AMA' filesep 'OUT' filesep];
        if ~Dir.exist(dire)
            mkdir(dire);
        end
        fname=[dire alias '_' dt '.mat'];
    end
    function fname=toFName(alias);
        fname=[Env.var('DATA') 'AMA' filesep alias '.mat'];
    end
    function alias=getAlias(stAlias,normType,fitType,bMeanCon,amaType)
        alias=[stAlias '_' normType '_' fitType '_' num2str(bMeanCon) '_' amaType];
    end
    function [obj,bSuccess]=load(aliasORfname,bRaw)
        if nargin < 2 || isempty(bRaw)
            bRaw=false;
        end
        if Fil.is(aliasORfname)
            fname=aliasORfname;
            [~,alias]=Fil.parts(fname);
        elseif contains(aliasORfname,filesep)
            error('File %s does not exist',aliasORfname);
        else
            alias=aliasORfname;
            alias=regexprep(alias,'\.mat$','');
            fname=AMA.alias2fname(alias);
            if ~Fil.exist(fname)
                if nargout > 1
                    obj=[];
                    bSuccess=0;
                    return
                else
                    error('File for alias %s does not exist',alias);
                end
            end
        end
        bSuccess=true;

        obj=AMA();
        load(fname);
        if bRaw
            obj=ama;
            return
        end
        flds=fieldnames(ama);
        for i = 1:length(flds)
            if isprop(obj,flds{i})
                obj.(flds{i})=ama.(flds{i});
            end
        end
        if nargout < 1
            assignin('base',['objAma_' alias],obj);
            assignin('base','obj',obj);
        end
        obj.loadTrainStim(); %should do automatically
        obj.loadedAlias=alias;
    end
    function opts=loadOpts(alias)
        name=sprintf('D_AMA_%s',alias);
        fname=which(name);
        if isempty(fname);
            fname=which([name '.cfg']);
        end
        if isempty(fname)
            error('AMA definition file with alias %s, does not exist',alias);
        end
        opts=Cfg.read(fname);
    end
end
methods
    function obj=AMALoader()
    end
end
end
