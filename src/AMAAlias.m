classdef AMAAlias < handle
properties
end
methods(Static)
end
methods
    function Alias(obj)
    end
    function Alias=gen(obj)
        NHash=obj.getNeuralHash;
        if ~isempty(NHash)
            NHash=['_' NHash];
        end
        if obj.rndSd==1
            SdStr='';
        else
            SdStr=['_' num2str(obj.randSd)];
        end
        if isempty(obj.bNoise) || ~obj.bNoise
            noiseStr='';
        else
            noiseStr='_Ns';
        end
        alias=[ ...
                  obj.stAlias '_' ...
                  obj.normType '_' ...
                  obj.fitType '_'  ...
                  num2str(obj.bMeanCon) '_' ...
                  obj.amaType ...
                  noiseStr ...
                  SdStr...
                  NHash;
        ];
        if nargout < 1
            if strcmp(obj.alias,alias)
                %error(['Alias already named ' alias '.']);
            elseif obj.bSet && bHasFit
                error('Fit data exists! Cannot change alias.')
            end
            obj.alias=alias;
        else
            Alias=alias;
        end
    end
    function hash=getNeuralHash(obj)
        [~,N]=AMA.getDN();
        N=rmfield(N,{'normType','bNoise'});
        flds=fieldnames(N);

        bHash=false;
        S=struct();
        for i = 1:size(flds,1)
            fld=flds{i,1};
            if ~isempty(obj.(fld)) && isfield(N,fld) && ~isequal(obj.(fld),N.(fld))
                S.(fld)=obj.(fld);
                bHash=true;
                fld
            end
        end
        if bHash
            hash=DataHash(S);
        else
            hash='';
        end
    end
end
end
