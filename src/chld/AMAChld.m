classdef AMAChld < handle
methods(Static)
end
properties(Hidden,Abstract)
    AMADeps
end
properties(Access=protected)
    AMA
    name
    prps

    fAMASet=0
    % hprps = things for copy() to copy
    % nprps = things to not copy()
end
methods(Static)
    function P=getP()
        P={...
        };
    end
end
methods
    function P=getPSave(obj)
        P=obj.getP();
    end
    function obj=AMAChld(varargin)
        if nargin > 1
            return
        end
    end
    function parse(obj)
        % for direct parsing
    end
    function init(obj)
    end
    function assert(obj)
    end
    function nbj=copy(obj,hprps,name)
        if nargin < 2
            hprps={};
        end
        if isprop(obj,'hprps')
            hprps=[obj.hprps hprps];
        end
        prps=obj.prps;
        if isprop(obj,'nprps')
            prps(ismember(prps,obj.nprps))=[];
        end
        if nargin >= 3 && ~isempty(name)
            name=name;
        else
            name=obj.name;
        end

        nbj=eval([name '();']);
        for i = 1:length(prps)
            prp=prps{i};
            nbj.(prp)=obj.(prp);
        end

        nbj.AMA=obj.AMA;

        if isempty(hprps)
            return
        end
        for i = 1:length(hprps)
            prp=prps{i};
            nbj.(prp)=obj.(prp);
        end

    end
end
methods(Hidden,Access=?AMA)
    function ama_receive(obj,AMA)
        obj.fAMASet=1;
        obj.AMA=AMA;

        % name
        obj.name=class(obj);

        % Get opts
        %if isfield(obj.AMA.Opts,obj.name)
        %    obj.Opts=obj.AMA.Opts.(name);
        %end

        % get own properties
        obj.prps=Obj.props(obj);

        % get other req objects
        %sisters=obj.AMA.ls_child;
        %for i = 1:length(sisters)
        %    sis=sisters{i};
        %    if any(strcmp(sis,obj.AMADeps))
        %        obj.(sis)=obj.AMA.(sis);
        %    end
        %end
    end
    function ama_apply_opts(obj,opts)
        Args.apply(obj,opts);
    end
end
end
